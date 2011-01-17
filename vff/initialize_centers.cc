#include "LaDaConfig.h"

#include <cmath>

#include <crystal/divide_and_conquer.h>
#include <crystal/neighbors.h>
#include <math/misc.h>

#include "vff.h"
  
namespace LaDa
{
  namespace vff
  { 
    void Vff :: first_neighbors_( std::vector< std::vector< math::rVector3d > >& _fn )
    {
      const size_t Nsites( structure.lattice->sites.size() );
      LADA_DO_NASSERT( Nsites != 2, "Expected two sites for VFF.\n" )
      typedef std::vector< std::vector< math::rVector3d > > t_FirstNeighbors;
      _fn.resize( structure.lattice->sites.size() );
      foreach( const Crystal::Lattice::t_Site &site, structure.lattice->sites )
        _fn[0].push_back( site.pos );
      _fn[1] = _fn[0];
      std::swap( _fn[1][0], _fn[1][1] ); 
      for( size_t i(0); i < Nsites; ++i )
        Crystal::find_first_neighbors( _fn[i], structure.lattice->cell, 4 );
    }

    bool Vff :: initialize_centers(bool _verbose)
    {
      LADA_ASSERT(structure.atoms.size() != 0, "Structure has no atoms.\n" )
      LADA_ASSERT(not math::is_zero(structure.cell.determinant()), "Structure with zero volume.\n")
      centers_.clear();
      { // Creates a list of centers
        t_Atoms :: iterator i_atom = structure.atoms.begin();
        t_Atoms :: iterator i_atom_end = structure.atoms.end();
        for(types::t_unsigned index=0; i_atom != i_atom_end; ++i_atom, ++index )
          centers_.push_back( AtomicCenter( structure, *i_atom, index ) );
      }

      // finds first neighbors on ideal lattice.
      t_FirstNeighbors fn;
      first_neighbors_( fn );

      // Tries to guess size of divide and conquer.
      const math::iVector3d nboxes( Crystal::guess_dnc_params( structure, 30 ) );
      types::t_real n(   structure.atoms.size()
                       / types::t_real( nboxes(0) * nboxes(1) * nboxes(2) ) );
       if( _verbose )
       {
         LADA_ROOT
         (
           comm,
           std::cout << "Will divide into " << nboxes(0) << "x"
                     << nboxes(1) << "x" << nboxes(2)
                     << " boxes of " << n << " atoms each.\n";
         )
       }
      // Then creates boxes.
      const types::t_real odist( 1.5e0 * std::sqrt( fn[0].front().squaredNorm() ) );
      Crystal::t_ConquerBoxes<t_Atom::t_Type> :: shared_ptr boxes
      (
        Crystal :: divide_and_conquer_boxes( structure, nboxes, odist )
      );
      // Finally calls algorithm.
      typedef Crystal::t_ConquerBoxes<t_Atom::t_Type>::type::const_iterator t_cit;
      t_cit i_box = boxes->begin();
      t_cit i_box_end = boxes->end();
      bool result( true );
      for(; result and i_box != i_box_end; ++i_box )
        result = build_tree_sort_dnc_( *i_box, fn );
      if( not result ) return false;
      

#     ifdef LADA_DEBUG
        check_tree(); 
#     endif
      if( _verbose ) 
      {
        LADA_ROOT( comm, std::cout << "First Neighbor tree successfully created.\n"; )
      }
      return true;
    } // Vff :: construct_bond_list

  } // namespace vff
} // namespace LaDa
