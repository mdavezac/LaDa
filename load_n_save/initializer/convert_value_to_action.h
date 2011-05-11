//
//  Version: $Id: convert_value_to_action.h 1226 2009-07-13 06:28:01Z davezac $
//

#ifndef _LADA_LOADNSAVE_INITIALIZER_CONVERT_VALUE_TO_ACTION_H_
#define _LADA_LOADNSAVE_INITIALIZER_CONVERT_VALUE_TO_ACTION_H_


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/conversion/cast.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/xpressive/regex_algorithms.hpp>
#include <boost/xpressive/regex_compiler.hpp>
#include <boost/mpl/bool.hpp>

#include <opt/fuzzy.h>
#include "convert_value_to_regex.h"
#include "parse_value.h"


namespace LaDa 
{
  namespace load_n_save
  {
    namespace initializer
    {
      //! Does nothing.
      inline bool convert_value_to_action( std::string const &_value,
                                           boost::mpl::bool_<false> const &) { return true; }
      //! Converts a string to an arithmetic action.
      template< class T_ACTION >
        typename boost::enable_if
        <
          typename boost::is_arithmetic<T_ACTION>::type, bool
        > :: type
        convert_value_to_action( std::string const &_value, T_ACTION &_action )
          { _action = boost::lexical_cast<T_ACTION>( _value ); return true; }
      //! Converts a value to a string
      inline bool convert_value_to_action( std::string const &_value, std::string &_action)
        { _action =_value; return true; }
      //! Converts a string to a string
      inline bool convert_value_to_action( std::string const &_value, const char *_action)
        { _action = _value.c_str(); return true; }
      //! Converts a string to a bool
      inline bool convert_value_to_action( std::string const &_value, bool &_action)
      {
        namespace bx = boost::xpressive;
        bx::sregex const re = bx::sregex::compile( convert_value_to_regex( true ) );
        _action = bx::regex_match( _value, re );
        return true;
      }

      template< class T_TYPE >
        typename boost::enable_if
        <
          typename is_action<T_TYPE>::type, bool
        >::type convert_value_to_action( std::string const &_value, T_TYPE& _action )
        {
          typename T_TYPE :: arg_type arg;
          if( not parse_value( _value, arg ) ) return false;
          return _action(arg);
        };

    }
  }
}

#endif
