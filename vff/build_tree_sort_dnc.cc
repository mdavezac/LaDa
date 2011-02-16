//
//  Version: $Id: build_tree_sort.cc 1079 2009-04-29 23:34:09Z Mayeul $
//
#include "LaDaConfig.h"

#include <algorithm>
#include <functional>

#include <boost/tuple/tuple.hpp>

#include <crystal/divide_and_conquer.h>

#include "exceptions.h"
#include "vff.h"
  
namespace LaDa
{
  namespace vff
  { 
    bool Vff :: build_tree_sort_dnc_( const Crystal::ConquerBox<t_Atom::t_Type>& _dnc, 
                                      const t_FirstNeighbors & _fn )
    {
      namespace bt = boost::tuples;
      typedef Crystal::ConquerBox<t_Atom::t_Type> :: t_States t_States;
      const t_States :: const_iterator i_state_begin = _dnc.begin();
      const t_States :: const_iterator i_state_end = _dnc.end();
      const math::rMatrix3d invcell(structure.cell.inverse());

      // loops over all centers in small box.
      for( t_States :: const_iterator i_center( i_state_begin );
           i_center != i_state_end; ++i_center )
      {
        if( not bt::get<1>( *i_center ) ) continue; // skip states outside small box.

        const size_t center_index( bt::get<0>( *i_center ) );
        t_Center &center( centers_[ center_index ] );
        const size_t site( center.atom().site );
        LADA_BASSERT( site < structure.lattice->sites.size(),
                      exceptions::site_index() << exceptions::integer(site) );
        const size_t neigh_site( site == 0 ? 1: 0 );
        const types::t_real cutoff = types::t_real(0.25) * _fn[site].front().squaredNorm();
                   
        // loops over all potential bonds (in large box)
        for( t_States :: const_iterator i_bond( i_state_begin );
             i_bond != i_state_end; ++i_bond )
        {
          const size_t bond_index( bt::get<0>( *i_bond ) );
          if( bond_index == center_index ) continue;

          const t_Center &bond( centers_[ bond_index ] );
          if( site == bond.atom().site ) continue;
          
          std::vector<math::rVector3d> :: const_iterator i_neigh = _fn[site].begin();
          const std::vector<math::rVector3d> :: const_iterator i_neigh_end = _fn[site].end();
          for(; i_neigh != i_neigh_end; ++i_neigh )
          {
            const math::rVector3d image
            ( 
              center.i_atom_->pos + *i_neigh - bond.i_atom_->pos
            );
            const math::rVector3d frac_image( invcell * image );
            const math::rVector3d frac_centered
            ( 
              frac_image[0] - rint( frac_image[0] ),
              frac_image[1] - rint( frac_image[1] ),
              frac_image[2] - rint( frac_image[2] )
            );
            const math::rVector3d cut( structure.cell * frac_centered );

            if( cut.squaredNorm() > cutoff ) continue;
            
            t_Centers :: iterator ii_bond( centers_.begin() + bond_index );
            center.bonds.push_back( t_Center ::__make__iterator__( ii_bond ) );
            const math::rVector3d trans
            (
              rint( frac_image[0] ),
              rint( frac_image[1] ),
              rint( frac_image[2] ) 
            );
            center.translations.push_back( trans );
            center.do_translates.push_back( not math::is_zero(trans.squaredNorm()) );

            if( center.bonds.size() == 4 ) break;
          } // loop over neighbors
          if( center.bonds.size() == 4 ) break;
        } // loop over bonds
      } // loop over centers_ 

      return true;
    } // Vff :: build_tree_sort_dnc

  } // namespace vff
} // namespace LaDa
