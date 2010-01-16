//
//  Version: $Id$
//
#ifndef LADA_MATH_IPOW_H
#define LADA_MATH_IPOW_H

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <opt/types.h>
#include "fuzzy.h"

namespace LaDa
{

  namespace math
  { 
    //! Exponentiation for integer powers.
    template<class T>
      inline T ipow(T base, size_t exponent) 
      {
        T result(1);
        for(; exponent; --exponent, result *= base);
        return result;
      }

    //! True if the number is an integer.
    inline bool is_integer(types::t_real const &x) { return is_zero(x - std::floor(x+0.1)); }
    //! True if the vector is composed of integers.
    inline bool is_integer(Eigen::Vector3d const &x)
      { return is_integer(x.x()) and is_integer(x.z()) and is_integer(x.y()); }
    //! True if the matrix is composed of integers.
    inline bool is_integer(Eigen::Matrix3d const &x)
    {
      for(size_t i(0); i < 3; ++i)
        for(size_t j(0); j < 3; ++j)
          if( not is_integer(x(i,j)) ) return false;
      return true;
    }

    //! True if two vectors are periodic images with respect to a cell.
    inline bool are_periodic_images( Eigen::Vector3d const &_a,
                                     Eigen::Vector3d const &_b, 
                                     Eigen::Matrix3d const &_inv_cell )
    {
      Eigen::Vector3d v = _inv_cell * (_a - _b); 
      return is_integer(v); 
    }

  } // namespace math

} // namespace LaDa

#endif
