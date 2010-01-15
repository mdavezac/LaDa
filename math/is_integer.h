//
//  Version: $Id$
//
#ifndef _INTEGER__H_
#define _INTEGER__H_

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include <cmath>

#include <opt/types.h>
#include <math/fuzzy.h>

namespace LaDa 
{ 

  namespace math
  {
    inline bool is_integer(types::t_real x) 
      { return Fuzzy::is_zero(x-std::floor(x+0.1e0)); }

    inline bool is_integer(Eigen::Vector3d const &x) 
    {
      for( size_t i(0); i < 3; ++i )
        if( not is_integer(x(i)) ) return false;
      return true;
    }

    inline bool is_integer(Eigen::Matrix3d const &x) 
    {
      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
          if( not is_integer(x(i, j)) ) return false;
      return true;
    }
  } // namespace math

} // namespace LaDa

#endif
