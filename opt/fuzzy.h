//
//  Version: $Id$
//
#ifndef _OPT_FUZZY_H_
#define _OPT_FUZZY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "types.h"

#ifdef _MPI
#include <mpi/mpi_object.h>
#endif

namespace opt 
{
  //! \brief Defines fuzzy math ordering operator for \a T_ARG
  //! \details if \a T_ARG is an integer type, signed or unsigned, then
  //! Fuzzy::less, Fuzzy::greater, Fuzzy::equal are exactly equivalent to <, >,
  //! and ==. If on the other hand \a T_ARG is types::real, then Fuzzy::less,
  //! Fuzzy::greater, Fuzzy::equal are defined with a fuzzy limit (types::tolerance)
  //! \param T_ARG is expected to be a scalar type
  template< class T_ARG>
    struct Fuzzy {
      //! \brief returns true if \f$ \_a  < \_b \f$
      //! \details if \a T_ARG is a types::real, return true if 
      //! \f$|\_a - \_b| < \f$ types::tolerance or\f$ \_a  < \_b \f$
      static bool leq( const T_ARG _a, const T_ARG _b ) 
        { return _a <= _b; }
      //! \brief returns true if \f$ \_a  > \_b \f$
      //! \details if \a T_ARG is a types::real, return true if 
      //! \f$|\_a - \_b| < \f$ types::tolerance or \f$ \_a  > \_b \f$
      static bool geq( const T_ARG _a, const T_ARG _b ) 
        { return _a >= _b; }
      //! \brief returns true if \f$ \_a  < \_b \f$
      //! \details if \a T_ARG is a types::real, return true if 
      //! \f$|\_a - \_b| > \f$ types::tolerance, and \f$ \_a  < \_b \f$
      static bool less( const T_ARG _a, const T_ARG _b ) 
        { return _a < _b; }
      //! \brief returns true if \f$ \_a  > \_b \f$
      //! \details if \a T_ARG is a types::real, return true if 
      //! \f$|\_a - \_b| > \f$ types::tolerance, and \f$ \_a  > \_b \f$
      static bool greater( const T_ARG _a, const T_ARG _b ) 
        { return _a > _b; }
      //! \brief returns true if \f$ \_a  == \_b \f$
      //! \details if \a T_ARG is a types::real, return true if 
      //! \f$|\_a - \_b| > \f$ types::tolerance, and \f$ \_a  > \_b \f$
      static bool equal( const T_ARG _a, const T_ARG _b ) 
        { return _a == _b; }
    };
  //! Specialize version for \a T_ARG a types::real, eg when fuzzyness is indeed needed 
  template<>
    struct Fuzzy<types::t_real> {
      //! \brief returns true if \f$ \_a  < \_b \f$
      //! \details if \a T_ARG is a types::real, return true if 
      //! \f$|\_a - \_b| < \f$ types::tolerance or\f$ \_a  < \_b \f$
      static bool leq( const types::t_real _a, const types::t_real _b ) 
        { return std::abs(_a - _b) < types::tolerance or _a < _b; }
      //! \brief returns true if \f$ \_a  > \_b \f$
      //! \details if \a T_ARG is a types::real, return true if 
      //! \f$|\_a - \_b| < \f$ types::tolerance or \f$ \_a  > \_b \f$
      static bool geq( const types::t_real _a, const types::t_real _b ) 
        { return std::abs(_a - _b) < types::tolerance or _a > _b; }
      //! \brief returns return true if 
      //! \f$|\_a - \_b| > \f$ types::tolerance, and \f$ \_a  < \_b \f$
      static bool less( const types::t_real _a, const types::t_real _b ) 
        { return std::abs(_a - _b) > types::tolerance and _a < _b; }
      //! \brief  return true if 
      //! \f$|\_a - \_b| > \f$ types::tolerance, and \f$ \_a  > \_b \f$
      static bool greater( const types::t_real _a, const types::t_real _b ) 
        { return std::abs(_a - _b) > types::tolerance and _a > _b; }
      //! \brief return true if \f$|\_a - \_b| > \f$ types::tolerance
      static bool equal( const types::t_real _a, const types::t_real _b ) 
        { return std::abs( _a - _b ) <  types::tolerance; }
    };


} // namspace opt


#endif
