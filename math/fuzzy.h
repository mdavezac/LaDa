//
//  Version: $Id$
//
#ifndef _OPT_FUZZY_H_
#define _OPT_FUZZY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>
#include "eigen.h"

namespace LaDa
{
  //! Namespace for mathematical primitives.
  namespace math
  {
    //! \brief returns true if \f$ \_a  < \_b \f$
    //! \details if \a T_ARG is a types::real, return true if 
    //! \f$|\_a - \_b| < \f$ types::tolerance or\f$ \_a  < \_b \f$
    template< class T_ARG >
    inline bool leq( const T_ARG _a, const T_ARG _b ) { return _a <= _b; }
    //! \brief returns true if \f$ \_a  > \_b \f$
    //! \details if \a T_ARG is a types::real, return true if 
    //! \f$|\_a - \_b| < \f$ types::tolerance or \f$ \_a  > \_b \f$
    template< class T_ARG >
    inline bool geq( const T_ARG _a, const T_ARG _b ) { return _a >= _b; }
    //! \brief returns true if \f$ \_a  < \_b \f$
    //! \details if \a T_ARG is a types::real, return true if 
    //! \f$|\_a - \_b| > \f$ types::tolerance, and \f$ \_a  < \_b \f$
    template< class T_ARG >
    inline bool le( const T_ARG _a, const T_ARG _b ) { return _a < _b; }
    //! \brief returns true if \f$ \_a  > \_b \f$
    //! \details if \a T_ARG is a types::real, return true if 
    //! \f$|\_a - \_b| > \f$ types::tolerance, and \f$ \_a  > \_b \f$
    template< class T_ARG >
    inline bool gt( const T_ARG _a, const T_ARG _b ) { return _a > _b; }
    //! \brief returns true if \f$ \_a  == \_b \f$
    //! \details if \a T_ARG is a types::real, return true if 
    //! \f$|\_a - \_b| > \f$ types::tolerance, and \f$ \_a  > \_b \f$
    template< class T_ARG >
    inline bool eq( const T_ARG _a, const T_ARG _b ) { return _a == _b; }
    //! \brief returns true if \f$ \_a  != \_b \f$
    //! \details if \a T_ARG is a types::real, return true if 
    //! \f$|\_a - \_b| >= \f$ types::tolerance.
    template< class T_ARG >
    inline bool neq( const T_ARG _a, const T_ARG _b ) { return _a != _b; }
    //! \brief returns true if \a _a  == 0.
    //! \details if \a T_ARG is a types::real, return true if 
    //!          \a _a < types::tolerance.
    template< class T_ARG > bool is_zero( T_ARG const &_a );

    //! Is a vector zero.
    template<> inline bool is_zero( rVector3d const &_a ) { return is_zero( _a.squaredNorm() );  }
    //! Is a vector zero.
    template<> inline bool is_zero( iVector3d const &_a )
    {
      if( _a(0) != 0 ) return false;
      if( _a(1) != 0 ) return false;
      return _a(2) == 0;
    }
    //! Is a matrix zero.
    template<> inline bool is_zero( rMatrix3d const &_a )
    {
      for(size_t i(0); i < 3; ++i)
        for(size_t j(0); j < 3; ++j)
          if( not is_zero(_a(i,j)) ) return false;
      return true;
    }
    template<class T> inline bool is_zero(T const &_a ) { return math::eq( _a, T(0) ); }

    //! \cond
    template<>
    inline bool leq<types::t_real>( const types::t_real _a,
                                    const types::t_real _b ) 
      { return std::abs(_a - _b) < types::tolerance or _a < _b; }
    template<>
    inline bool geq<types::t_real>( const types::t_real _a,
                                    const types::t_real _b ) 
      { return std::abs(_a - _b) < types::tolerance or _a > _b; }
    template<>
    inline bool le<types::t_real>( const types::t_real _a,
                                   const types::t_real _b ) 
      { return std::abs(_a - _b) > types::tolerance and _a < _b; }
    template<>
    inline bool gt<types::t_real>( const types::t_real _a,
                                   const types::t_real _b ) 
      { return std::abs(_a - _b) > types::tolerance and _a > _b; }
    template<>
    inline bool eq<types::t_real>( const types::t_real _a,
                                   const types::t_real _b ) 
      { return std::abs( _a - _b ) <  types::tolerance; }
    template<>
    inline bool neq<types::t_real>( const types::t_real _a,
                                    const types::t_real _b ) 
      { return std::abs( _a - _b ) >=  types::tolerance; }
    //! \endcond
  } // namepace math
} // namespace LaDa

#endif
