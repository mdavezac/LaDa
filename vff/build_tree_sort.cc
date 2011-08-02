#include "LaDaConfig.h"

#include <algorithm>
#include <functional>

#include "vff.h"
  
namespace LaDa
{
  namespace vff
  { 
    bool Vff :: build_tree_sort_(const t_FirstNeighbors & _fn)
    {
      t_Centers :: iterator i_begin = centers_.begin();
      t_Centers :: iterator i_end = centers_.end();
      t_Centers :: iterator i_center, i_bond;
      for( i_center = i_begin; i_center != i_end; ++i_center )
      {
        const size_t site( i_center->atom().site );
        LADA_DO_NASSERT( site > structure.lattice->sites.size(), "Unindexed site.\n" )
        const size_t neigh_site( site == 0 ? 1: 0 );
        const types::t_real cutoff = types::t_real(0.25) * _fn[site].front().squaredNorm();
                   

        for( i_bond = i_begin; i_bond != i_end; ++i_bond)
        {
          if( i_bond == i_center ) continue;
          if( i_bond->atom().site == site ) continue;
          
          std::vector<math::rVector3d> :: const_iterator i_neigh = _fn[site].begin();
          const std::vector<math::rVector3d> :: const_iterator i_neigh_end = _fn[site].end();
          for(; i_neigh != i_neigh_end; ++i_neigh )
          {
            const math::rVector3d image
            ( 
              i_center->i_atom_->pos + *i_neigh - i_bond->i_atom_->pos 
            );
            const math::rVector3d frac_image( (!structure.cell) * image );
            const math::rVector3d frac_centered
            ( 
              frac_image[0] - rint( frac_image[0] ),
              frac_image[1] - rint( frac_image[1] ),
              frac_image[2] - rint( frac_image[2] )
            );
            const math::rVector3d cut( structure.cell * frac_centered );

            if( cut.squaredNorm() > cutoff ) continue;
            
            i_center->bonds.push_back( t_Center ::__make__iterator__(i_bond) );
            const math::rVector3d trans
            (
              rint( frac_image[0] ),
              rint( frac_image[1] ),
              rint( frac_image[2] ) 
            );
            i_center->translations.push_back( trans );
            i_center->do_translates.push_back( not math::is_null(trans.squaredNorm()) );

            if( i_center->bonds.size() == 4 ) break;
          } // loop over neighbors
          if( i_center->bonds.size() == 4 ) break;
        } // loop over bonds
      } // loop over centers 

      return true;
    } // Vff :: build_tree_sort

  } // namespace vff
} // namespace LaDa
