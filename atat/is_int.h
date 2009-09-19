//
//  Version: $Id$
//

#ifndef LADA_ATAT_IS_INT_H
#define LADA_ATAT_IS_INT_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<cmath>
#include "vectmac.h"

namespace LaDa 
{
  namespace atat
  {
    // returns true if _a is an integer matrix.
    inline bool is_integer( atat::rMatrix3d const &_a, types::t_real _tolerance = 1e-12 )
    {
      for(size_t i(0); i < 3; ++i)
        for(size_t j(0); j < 3; ++j)
          if( std::abs( _a(i,j) - std::floor(_a(i,j)+0.1) ) < _tolerance ) return false;
      return true;
    }

    // returns true if _a is an integer vector.
    inline bool is_integer( atat::rVector3d const &_a, types::t_real _tolerance = 1e-12 )
    {
      for(size_t i(0); i < 3; ++i)
        if( std::abs( _a(i) - std::floor(_a(i)+0.1) ) < _tolerance ) return false;
      return true;
    }

  } // namespace atat

} // namespace LaDa



#endif
