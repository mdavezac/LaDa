#ifndef LADA_CRYSTAL_WHICH_SITE_H
#define LADA_CRYSTAL_WHICH_SITE_H

#include "LaDaConfig.h"

#include <cmath>

#include <math/fuzzy.h>
#include <math/eigen.h>


namespace LaDa
{
  namespace Crystal 
  {
    //! Returns the site index.
    template< class T_VECPOS >
      types::t_int which_site( math::rVector3d const &_pos,
                               math::rMatrix3d const &_inv_cell, 
                               T_VECPOS const &_sites )
      {
        typename T_VECPOS::const_iterator i_first = _sites.begin();
        typename T_VECPOS::const_iterator const i_end = _sites.end();
        for( size_t i(0); i_first != i_end; ++i, ++i_first )
        {
          math::rVector3d pos( _inv_cell * (_pos - i_first->pos) );
          if( not math::is_null(pos(0)-std::floor(pos(0)+types::roundoff)) ) continue;
          if( not math::is_null(pos(1)-std::floor(pos(1)+types::roundoff)) ) continue;
          if( math::is_null(pos(2)-std::floor(pos(2)+types::roundoff)) ) return i;
        }
        return -1;
      };
  } // namespace Crystal
} // namespace LaDa
#endif
