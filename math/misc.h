#ifndef LADA_MATH_IPOW_H
#define LADA_MATH_IPOW_H

#include "LaDaConfig.h"

#include <opt/types.h>
#include "eigen.h"
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
    inline bool is_integer(math::rVector3d const &x)
      { return is_integer(x.x()) and is_integer(x.z()) and is_integer(x.y()); }
    //! True if the matrix is composed of integers.
    inline bool is_integer(math::rMatrix3d const &x)
    {
      for(size_t i(0); i < 3; ++i)
        for(size_t j(0); j < 3; ++j)
          if( not is_integer(x(i,j)) ) return false;
      return true;
    }
    //! True if the obect is composed of integers.
    template< class T1, class T2, int T3>
      inline bool is_integer(Eigen::Product<T1, T2, T3> const &x)
        { return is_integer( x.eval() ); }


    //! True if two vectors are periodic images with respect to a cell.
    inline bool are_periodic_images( math::rVector3d const &_a,
                                     math::rVector3d const &_b, 
                                     math::rMatrix3d const &_inv_cell )
    {
      math::rVector3d v = _inv_cell * (_a - _b); 
      return is_integer(v); 
    }


    //! Defines a floor function.
    inline rMatrix3d floor(rMatrix3d const &_matrix)
    {
      rMatrix3d result;
      for(size_t i(0); i < 3; ++i)
        for(size_t j(0); j < 3; ++j)
          result(i,j) = std::floor(_matrix(i,j));
      return result;
    }
    //! Defines a floor function.
    inline math::rVector3d floor(math::rVector3d const &_v) 
    {
      return math::rVector3d( std::floor(_v(0)), 
                              std::floor(_v(1)), 
                              std::floor(_v(2)) );
    }

    //! Rounding off.
    inline types::t_real round(types::t_real r) 
      { return (r > 0.0) ? std::floor(r + 0.5) : std::ceil(r - 0.5); }

    //! Rounding off.
    inline math::rVector3d round(math::rVector3d const &_v)
      { return math::rVector3d( round(_v(0)), round(_v(1)), round(_v(2)) ); }

    //! Maximum component.
    template<class T> inline T max( Eigen::Matrix<T, 3, 1>const &_v)
      { return std::max(std::max(_v(0), _v(1)), _v(2)); }

  } // namespace math

} // namespace LaDa

#endif
