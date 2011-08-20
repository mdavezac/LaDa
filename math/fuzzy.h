#ifndef LADA_MATH_FUZZY_H
#define LADA_MATH_FUZZY_H

#include "LaDaConfig.h"

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/type_traits/is_floating_point.hpp>

#include <opt/types.h>
#include "eigen.h"

namespace LaDa
{
  //! Namespace for mathematical primitives.
  namespace math
  {
#   ifdef LADA_INTEGRAL
#     error LADA_INTEGRAL already defined.
#   endif
#   ifdef LADA_REAL
#     error LADA_REAL already defined.
#   endif
#   ifdef LADA_ARITH
#     error LADA_ARITH already defined.
#   endif
#   define LADA_INTEGRAL typename boost::enable_if<boost::is_integral<T_ARG>, bool> :: type
#   define LADA_REAL typename boost::enable_if<boost::is_floating_point<T_ARG>, bool> :: type
#   define LADA_ARITH typename boost::enable_if<boost::is_arithmetic<T_ARG>, bool> :: type
    //! \brief True if \f$|\_a - \_b| < \_tol\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_REAL eq(T_ARG const& _a, T_ARG const& _b, T_ARG const _tol = types::tolerance )
      { return std::abs(_a-_b) < _tol; }
    //! \brief True if \f$|\_a - \_b| < \_tol\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_INTEGRAL eq(T_ARG const& _a, T_ARG const& _b, T_ARG const& _tol)
      { return std::abs(_a-_b) < _tol; }
    //! \brief True if \f$\_a == \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_INTEGRAL eq(T_ARG const& _a, T_ARG const& _b) { return _a == _b; }

    //! \brief True if \f$|\_a - \_b| < \_tol\f$ or \f$\_a < \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_REAL leq(T_ARG const& _a, T_ARG const& _b, T_ARG const _tol = types::tolerance )
      { return _a <= _b or eq(_a, _b); }
    //! \brief True if \f$|\_a - \_b| < \_tol\f$ or \f$\_a < \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_INTEGRAL leq(T_ARG const& _a, T_ARG const& _b, T_ARG const& _tol)
      { return _a <= _b or eq(_a, _b); }
    //! \brief True if \f$\_a <= \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_INTEGRAL leq(T_ARG const& _a, T_ARG const& _b) { return _a <= _b; }

    //! \brief True if \f$|\_a - \_b| < \_tol\f$ or \f$\_a > \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_REAL geq(T_ARG const& _a, T_ARG const& _b, T_ARG const _tol = types::tolerance )
      { return _a >= _b or eq(_a, _b); }
    //! \brief True if \f$|\_a - \_b| < \_tol\f$ or \f$\_a > \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_INTEGRAL geq(T_ARG const& _a, T_ARG const& _b, T_ARG const& _tol)
      { return _a >= _b or eq(_a, _b); }
    //! \brief True if \f$\_a >= \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_INTEGRAL geq(T_ARG const& _a, T_ARG const& _b) { return _a >= _b; }

    //! \brief True if \f$|\_a - \_b| > \_tol\f$ and \f$\_a < \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_ARITH lt(T_ARG const& _a, T_ARG const& _b)   { return not geq(_a, _b); }
    //! \brief True if \f$|\_a - \_b| > \_tol\f$ and \f$\_a < \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_ARITH lt(T_ARG const& _a, T_ARG const& _b, T_ARG const& _tol)   { return not geq(_a, _b, _tol); }

    //! \brief True if \f$|\_a - \_b| > \_tol\f$ and \f$\_a > \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_ARITH gt(T_ARG const& _a, T_ARG const& _b)   { return not leq(_a, _b); }
    //! \brief True if \f$|\_a - \_b| > \_tol\f$ and \f$\_a > \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_ARITH gt(T_ARG const& _a, T_ARG const& _b, T_ARG const& _tol)   { return not leq(_a, _b, _tol); }

    //! \brief True if \f$|\_a - \_b| > \_tol\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_ARITH neq(T_ARG const& _a, T_ARG const& _b)   { return not eq(_a, _b); }
    //! \brief True if \f$|\_a - \_b| > \_tol\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_ARITH neq(T_ARG const& _a, T_ARG const& _b, T_ARG const& _tol)   { return not eq(_a, _b, _tol); }
    
    //! True if the number is an integer.
    template<class T_ARG>
    inline LADA_REAL is_integer(T_ARG const& x, T_ARG const _tol = types::tolerance)
        { return eq(x, std::floor(x+0.1), _tol); }
    //! True if the number is an integer.
    template<class T_ARG> inline LADA_INTEGRAL is_integer(T_ARG const&, T_ARG const&) { return true; }
    //! True if the number is an integer.
    template<class T_ARG> inline LADA_INTEGRAL is_integer(T_ARG const&) { return true; }

    //! \brief returns true if \a _a  == 0.
    //! \details if \a T_ARG is a types::real, return true if 
    //!          \a _a < types::tolerance.
    template<class T_ARG> inline LADA_ARITH is_null(T_ARG const& _a, T_ARG const& _tol)
      { return eq(_a, T_ARG(0), _tol); }
    //! \brief returns true if \a _a  == 0.
    //! \details if \a T_ARG is a types::real, return true if 
    //!          \a _a < types::tolerance.
    template<class T_ARG> inline LADA_ARITH is_null(T_ARG const& _a) { return eq(_a, T_ARG(0)); }

    //! \brief returns true if \a _a  == 0.
    //! \details if \a T_ARG is a types::real, return true if 
    //!          \a _a < types::tolerance.
    template<class T_ARG> inline LADA_ARITH is_identity(T_ARG const& _a, T_ARG const& _tol)
      { return eq(_a, T_ARG(1), _tol); }
    //! \brief returns true if \a _a  == 0.
    //! \details if \a T_ARG is a types::real, return true if 
    //!          \a _a < types::tolerance.
    template<class T_ARG> inline LADA_ARITH is_identity(T_ARG const& _a) { return eq(_a, T_ARG(1)); }
#   undef LADA_INTEGRAL
#   undef LADA_REAL
#   undef LADA_ARITH
  } // namepace math
} // namespace LaDa

#endif
