//
//  Version: $Id$
//
#ifndef _OPT_FUZZY_H_
#define _OPT_FUZZY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "types.h"

#include <mpi/mpi_object.h>



//! \brief Defines fuzzy math ordering operator for \a T_ARG
//! \details if \a T_ARG is an integer type, signed or unsigned, then
//!          Fuzzy::le, Fuzzy::leq, Fuzzy::ge, Fuzzy::geq, Fuzzy::eq are
//!          exactly equivalent to <, <=, >, =>, and ==. If on the other hand
//!          \a T_ARG is types::real, then Fuzzy::less, Fuzzy::greater,
//!          Fuzzy::equal are defined with a fuzzy limit (types::tolerance)
namespace Fuzzy
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
    template< class T_ARG >
    inline bool is_zero( const T_ARG _a ) { return Fuzzy::eq( _a, T_ARG(0) ); }
} // end namsepace Fuzzy


//! \cond
namespace Fuzzy
{
    template<>
    inline bool leq<types::t_real>( const types::t_real _a,
                                     types::t_real _b ) 
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
}
//! \endcond

#endif
