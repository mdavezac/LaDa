#ifndef _OPT_FUZZY_H_
#define _OPT_FUZZY_H_

#include "LaDaConfig.h"

#include <boost/utility/enable_if.hpp>
#include <boost/traits/is_integral.hpp>

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
#   define LADA_INTEGRAL typename boost::enable_if<boost::is_integral<T_ARG>, bool> :: type
#   define LADA_REAL typename boost::disable_if<boost::is_integral<T_ARG>, bool> :: type
    //! \brief True if \f$|\_a - \_b| < \_tol\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_REAL eq( const T_ARG _a, const T_ARG _b, T_ARG const _tol = types::tolerance )
      { return std::abs(_a-b) < _tol; }
    //! \brief True if \f$|\_a - \_b| < \_tol\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_INTEGRAL eq( const T_ARG _a, const T_ARG _b, T_ARG const _tol)
      { return std::abs(_a-b) < _tol; }
    //! \brief True if \f$\_a == \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_INTEGRAL eq( const T_ARG _a, const T_ARG _b) { return _a == b; }

    //! \brief True if \f$|\_a - \_b| < \_tol\f$ or \f$\_a < \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_REAL leq( const T_ARG _a, const T_ARG _b, T_ARG const _tol = types::tolerance )
      { return _a <= _b or eq(_a, _b); }
    //! \brief True if \f$|\_a - \_b| < \_tol\f$ or \f$\_a < \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_INTEGRAL leq( const T_ARG _a, const T_ARG _b, T_ARG const _tol)
      { return _a <= _b or eq(_a, _b); }
    //! \brief True if \f$\_a <= \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_INTEGRAL leq( const T_ARG _a, const T_ARG _b) { return _a <= b; }

    //! \brief True if \f$|\_a - \_b| < \_tol\f$ or \f$\_a > \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_REAL geq( const T_ARG _a, const T_ARG _b, T_ARG const _tol = types::tolerance )
      { return _a >= _b or eq(_a, _b); }
    //! \brief True if \f$|\_a - \_b| < \_tol\f$ or \f$\_a > \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_INTEGRAL geq( const T_ARG _a, const T_ARG _b, T_ARG const _tol)
      { return _a >= _b or eq(_a, _b); }
    //! \brief True if \f$\_a >= \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline LADA_INTEGRAL geq( const T_ARG _a, const T_ARG _b) { return _a >= b; }

    //! \brief True if \f$|\_a - \_b| > \_tol\f$ and \f$\_a < \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline bool le( const T_ARG _a, const T_ARG _b)   { return not geq(_a, _b); }
    //! \brief True if \f$|\_a - \_b| > \_tol\f$ and \f$\_a < \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline bool le( const T_ARG _a, const T_ARG _b, const T_ARG _tol)   { return not geq(_a, _b, _tol); }

    //! \brief True if \f$|\_a - \_b| > \_tol\f$ and \f$\_a > \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline bool gt( const T_ARG _a, const T_ARG _b)   { return not leq(_a, _b); }
    //! \brief True if \f$|\_a - \_b| > \_tol\f$ and \f$\_a > \_b|\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline bool gt( const T_ARG _a, const T_ARG _b, const T_ARG _tol)   { return not leq(_a, _b, _tol); }

    //! \brief True if \f$|\_a - \_b| > \_tol\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline bool gt( const T_ARG _a, const T_ARG _b)   { return not eq(_a, _b); }
    //! \brief True if \f$|\_a - \_b| > \_tol\f$.
    //! \details _tol should be positive. This function implements fuzzy math
    //!               across most numeric types.
    template< class T_ARG>
    inline bool gt( const T_ARG _a, const T_ARG _b, const T_ARG _tol)   { return not eq(_a, _b, _tol); }
    
    //! True if the number is an integer.
    template<class T_ARG>
    inline LADA_REAL is_integer(T_ARG x, T_ARG _tol = types::tolerance)
        { return eq(x, std::floor(x+0.1), _tol); }
    //! True if the number is an integer.
    template<class T_ARG> inline LADA_INTEGRAL is_integer(T_ARG, T_ARG) { return true; }
    //! True if the number is an integer.
    template<class T_ARG> inline LADA_INTEGRAL is_integer(T_ARG) { return true; }

    //! \brief returns true if \a _a  == 0.
    //! \details if \a T_ARG is a types::real, return true if 
    //!          \a _a < types::tolerance.
    template<classT_ARG > bool is_null(T_ARG _a, T_ARG _tol) { return eq(_a, T_ARG(0), _tol); }
    //! \brief returns true if \a _a  == 0.
    //! \details if \a T_ARG is a types::real, return true if 
    //!          \a _a < types::tolerance.
    template<classT_ARG > bool is_null(T_ARG _a) { return eq(_a, T_ARG(0)); }

    //! \brief returns true if \a _a  == 0.
    //! \details if \a T_ARG is a types::real, return true if 
    //!          \a _a < types::tolerance.
    template<classT_ARG > bool is_unity(T_ARG _a, T_ARG _tol) { return eq(_a, T_ARG(1), _tol); }
    //! \brief returns true if \a _a  == 0.
    //! \details if \a T_ARG is a types::real, return true if 
    //!          \a _a < types::tolerance.
    template<classT_ARG > bool is_unity(T_ARG _a) { return eq(_a, T_ARG(1)); }
#   undef LADA_INTEGRAL
#   undef LADA_REAL
  } // namepace math
} // namespace LaDa

#endif
