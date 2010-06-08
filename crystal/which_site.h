//
//  Version: $Id$
//
#ifndef LADA_CRYSTAL_WHICH_SITE_H
#define LADA_CRYSTAL_WHICH_SITE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cmath>

#include <math/fuzzy.h>
#include <math/eigen.h>
#include <limits>


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
        const types::t_real roundoff = 1e1 * std::numeric_limits<types::t_real>::epsilon();
        typename T_VECPOS::const_iterator i_first = _sites.begin();
        typename T_VECPOS::const_iterator const i_end = _sites.end();
        for( size_t i(0); i_first != i_end; ++i, ++i_first )
        {
          math::rVector3d pos( _inv_cell * (_pos - i_first->pos) );
          if( not math::is_zero(pos(0)-std::floor(pos(0)+roundoff)) ) continue;
          if( not math::is_zero(pos(1)-std::floor(pos(1)+roundoff)) ) continue;
          if( math::is_zero(pos(2)-std::floor(pos(2)+roundoff)) ) return i;
          if( math::is_zero(pos(2)) ) return i;
        }
        return -1;
      };
  } // namespace Crystal
} // namespace LaDa
#endif
