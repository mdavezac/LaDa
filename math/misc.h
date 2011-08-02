#ifndef LADA_MATH_IPOW_H
#define LADA_MATH_IPOW_H

#include "LaDaConfig.h"

#include <boost/numeric/conversion/converter.hpp>
#include <limits>

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

#   ifdef LADA_MACRO
#     error LADA_MACRO already defined.
#   endif
#   define LADA_MACRO(name, type, default_) \
      inline bool name(type const &x, default_ ) \
        { return name(x(0), _tol) and name(x(1), _tol) and name(x(2), _tol); }   \
      inline bool name(math::rVector3d const &x) \
        { return name(x(0)) and name(x(1)) and name(x(2)); }   
    LADA_MACRO(is_integer,  math::rVector3d, types::t_real _tol = types::tolerance);
    LADA_MACRO(is_null,     math::rVector3d, types::t_real _tol = types::tolerance);
    LADA_MACRO(is_identity, math::rVector3d, types::t_real _tol = types::tolerance);
    LADA_MACRO(is_integer,  math::iVector3d, types::t_int _tol = 0);
    LADA_MACRO(is_null,     math::iVector3d, types::t_int _tol = 0);
    LADA_MACRO(is_identity, math::iVector3d, types::t_int _tol = 0);
#   undef LADA_MACRO
#   define LADA_MACRO(name, type, default_) \
      inline bool name(type const &x, default_)         \
      {                                                 \
        for(size_t i(0); i < 3; ++i)                    \
          for(size_t j(0); j < 3; ++j)                  \
            if( not name(x(i,j), _tol) ) return false;  \
        return true;                                    \
      }                                                 \
      inline bool name(type const &x)                   \
      {                                                 \
        for(size_t i(0); i < 3; ++i)                    \
          for(size_t j(0); j < 3; ++j)                  \
            if( not name(x(i,j)) ) return false;        \
        return true;                                    \
      }
    LADA_MACRO(is_integer, math::rMatrix3d, types::t_real _tol = types::tolerance)
    LADA_MACRO(is_null,    math::rMatrix3d, types::t_real _tol = types::tolerance)
    LADA_MACRO(is_integer, math::iMatrix3d, types::t_int _tol = 0)
    LADA_MACRO(is_null,    math::iMatrix3d, types::t_int _tol = 0)
#   undef LADA_MACRO
#   define LADA_MACRO(name, type, default_) \
      inline bool name(type const &x, default_)         \
      {                                                 \
        for(size_t i(0); i < 3; ++i)                    \
          for(size_t j(0); j < 3; ++j)                  \
            if( i ==j and not name(x(i,j), _tol) ) return false;  \
           else if( not is_null(x(i,j), _tol) ) return false;  \
        return true;                                    \
      }                                                 \
      inline bool name(type const &x)                   \
      {                                                 \
        for(size_t i(0); i < 3; ++i)                    \
          for(size_t j(0); j < 3; ++j)                  \
            if( i == j not name(x(i,j)) ) return false; \
            if( not is_null(x(i,j)) ) return false;     \
        return true;                                    \
      }
    LADA_MACRO(is_unity, math::rMatrix3d, types::t_real _tol = types::tolerance)
    LADA_MACRO(is_unity, math::iMatrix3d, types::t_int _tol = 0)
#   undef LADA_MACRO

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

    //! Casts to lower integer accounting for numerical noise.
    template<int D> inline Eigen::Matrix<types::t_int, D, 1> 
      floor_int( Eigen::Matrix<types::t_real, D, 1> const &_t )
      { 
        typedef boost::numeric::converter<types::t_int,types::t_real> converter;
        Eigen::Matrix<types::t_int, D, 1> result;
        for(int i(0); i < D; ++i)
          result(i) = converter::nearbyint( std::floor(_t(i) + types::roundoff) );
             
        return result;
      }

  } // namespace math

} // namespace LaDa

#endif
