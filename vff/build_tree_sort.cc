//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <algorithm>
#include <functional>

#include <opt/ndim_iterator.h>
#include <opt/atat.h>

#include "vff.h"
  
namespace LaDa
{
  namespace Vff
  { 
    bool Vff :: build_tree_sort_(const t_FirstNeighbors & _fn)
    {
      t_Centers :: iterator i_begin = centers.begin();
      t_Centers :: iterator i_end = centers.end();
      t_Centers :: iterator i_center, i_bond;
      for( i_center = i_begin; i_center != i_end; ++i_center )
      {
        const size_t site( i_center->Origin().site );
        __DOASSERT( site > structure.lattice->sites.size(), "Unindexed site.\n" )
        const size_t neigh_site( site == 0 ? 1: 0 );
        const types::t_real cutoff = types::t_real(0.25) * atat::norm2( _fn[site].front() );
                   

        for( i_bond = i_begin; i_bond != i_end; ++i_bond)
        {
          if( i_bond == i_center ) continue;
          if( i_bond->Origin().site == site ) continue;
          
          std::vector<atat::rVector3d> :: const_iterator i_neigh = _fn[site].begin();
          const std::vector<atat::rVector3d> :: const_iterator i_neigh_end = _fn[site].end();
          for(; i_neigh != i_neigh_end; ++i_neigh )
          {
            const atat::rVector3d image
            ( 
              i_center->origin->pos + *i_neigh - i_bond->origin->pos 
            );
            const atat::rVector3d frac_image( (!structure.cell) * image );
            const atat::rVector3d frac_centered
            ( 
              frac_image[0] - rint( frac_image[0] ),
              frac_image[1] - rint( frac_image[1] ),
              frac_image[2] - rint( frac_image[2] )
            );
            const atat::rVector3d cut( structure.cell * frac_centered );

            if( atat::norm2( cut ) > cutoff ) continue;
            
            i_center->bonds.push_back( t_Center ::__make__iterator__(i_bond) );
            const atat::rVector3d trans
            (
              rint( frac_image[0] ),
              rint( frac_image[1] ),
              rint( frac_image[2] ) 
            );
            i_center->translations.push_back( trans );
            i_center->do_translates.push_back
            ( 
              atat::norm2(trans) > atat::zero_tolerance 
            );

            if( i_center->bonds.size() == 4 ) break;
          } // loop over neighbors
          if( i_center->bonds.size() == 4 ) break;
        } // loop over bonds
      } // loop over centers 

      return true;
    } // Vff :: build_tree_sort

  } // namespace vff
} // namespace LaDa
