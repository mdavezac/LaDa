//
//  Version: $Id$
//
#ifndef LADA_CRYSTAL_WHICH_SITE_H
#define LADA_CRYSTAL_WHICH_SITE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cmath>

#include <atat/vectmac.h>
#include <atat/machdep.h>
#include <opt/fuzzy.h>


namespace LaDa
{
  namespace Crystal 
  {
    //! Returns the site index.
    template< class T_VECPOS >
      types::t_int which_site( atat::rVector3d const &_pos,
                               atat::rMatrix3d const &_inv_cell, 
                               T_VECPOS const &_sites )
      {
        typename T_VECPOS::const_iterator i_first = _sites.begin();
        typename T_VECPOS::const_iterator const i_end = _sites.end();
        for( size_t i(0); i_first != i_end; ++i, ++i_first )
        {
          atat::rVector3d pos( _inv_cell * (_pos - i_first->pos) );
          pos(0) -= std::floor(pos(0)+0.01);
          if( not Fuzzy::is_zero(pos(0)) ) continue;
          pos(1) -= std::floor(pos(1)+0.01);
          if( not Fuzzy::is_zero(pos(1)) ) continue;
          pos(2) -= std::floor(pos(2)+0.01);
          if( Fuzzy::is_zero(pos(2)) ) return i;
        }
        return -1;
      };
  } // namespace Crystal
} // namespace LaDa
#endif
