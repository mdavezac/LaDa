//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <cmath>

#include <crystal/ideal_lattice.h>
#include <crystal/divide_and_conquer.h>
#include <math/misc.h>

#include "vff.h"
  
namespace LaDa
{
  namespace Vff
  { 
    void Vff :: first_neighbors_( std::vector< std::vector< math::rVector3d > >& _fn )
    {
      const size_t Nsites( structure.lattice->sites.size() );
      __DOASSERT( Nsites != 2, "Expected two sites for VFF.\n" )
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
      centers.clear();
      { // Creates a list of centers
        t_Atoms :: iterator i_atom = structure.atoms.begin();
        t_Atoms :: iterator i_atom_end = structure.atoms.end();
        for(types::t_unsigned index=0; i_atom != i_atom_end; ++i_atom, ++index )
          centers.push_back( AtomicCenter( structure, *i_atom, index ) );
      }

      // finds first neighbors on ideal lattice.
      t_FirstNeighbors fn;
      first_neighbors_( fn );

//     // Checks if ideal structure. In which case use smith normal tree building.
//     if( math::is_integer(structure.lattice->cell.inverse() * structure.cell) )
//     {
//       if( _verbose ) 
//       { 
//         __ROOTCODE
//         ( 
//           MPI_COMM,
//           std::cout << "Trying to create first neighbor tree "
//                        "using smith normal form algorithm.\n";
//         )
//       }
//       if( build_tree_smith_( fn ) )
//       {
//         __DODEBUGCODE( check_tree(); )
//         if( _verbose )
//         { 
//           __ROOTCODE( MPI_COMM,
//                       std::cout << "First Neighbor tree successfully created.\n"; ) 
//         }
//         return true;
//       }
//       __ROOTCODE( MPI_COMM, std::cout << "Failed.\n"; )
//       t_Centers :: iterator i_center = centers.begin();
//       t_Centers :: iterator i_center_end = centers.end();
//       for(; i_center != i_center_end; ++i_center ) i_center->bonds.clear();
//     } 
//
//     const size_t Nperbox( 10 );
//     if( structure.atoms.size() < Nperbox )
//     {
//       if( _verbose )
//       {
//         __ROOTCODE
//         ( 
//           MPI_COMM,
//           std::cout << "Creating first neighbor tree using standard algorithm.\n";
//         )
//       }
//       if( not build_tree_sort_( fn ) ) return false;
//       __DODEBUGCODE( check_tree(); )
//       if( _verbose )
//       {
//         __ROOTCODE( MPI_COMM,
//                     std::cout << "First Neighbor tree successfully created.\n"; )
//       }
//       return true;
//     }
//      
//     if( _verbose )
//     {
//       __ROOTCODE
//       (
//         MPI_COMM,
//         std::cout << "Creating first neighbor tree using "
//                      "divide-and-conquer algorithm.\n";
//       )
//     }
      // Tries to guess size of divide and conquer.
      const math::iVector3d nboxes( Crystal::guess_dnc_params( structure, 30 ) );
      types::t_real n(   structure.atoms.size()
                       / types::t_real( nboxes(0) * nboxes(1) * nboxes(2) ) );
       if( _verbose )
       {
        __ROOTCODE
        (
          MPI_COMM,
          std::cout << "Will divide into " << nboxes(0) << "x"
                    << nboxes(1) << "x" << nboxes(2)
                    << " boxes of " << n << " atoms each.\n";
        )
       }
      // Then creates boxes.
      const types::t_real odist( 1.5e0 * std::sqrt( fn[0].front().squaredNorm() ) );
      Crystal::t_ConquerBoxes<types::t_real> :: shared_ptr boxes
      (
        Crystal :: divide_and_conquer_boxes( structure, nboxes, odist )
      );
      // Finally calls algorithm.
      typedef Crystal::t_ConquerBoxes<types::t_real>::type::const_iterator t_cit;
      t_cit i_box = boxes->begin();
      t_cit i_box_end = boxes->end();
      bool result( true );
      for(; result and i_box != i_box_end; ++i_box )
        result = build_tree_sort_dnc_( *i_box, fn );
      if( not result ) return false;
      

      __DODEBUGCODE( check_tree(); )
      if( _verbose ) 
      {
        __ROOTCODE( MPI_COMM, std::cout << "First Neighbor tree successfully created.\n"; )
      }
      return true;
    } // Vff :: construct_bond_list

  } // namespace vff
} // namespace LaDa
