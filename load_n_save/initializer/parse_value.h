//
//  Version: $Id: parse_value.h 1226 2009-07-13 06:28:01Z davezac $
//

#ifndef _LADA_LOADNSAVE_INITIALIZER_PARSE_VALUE_H_
#define _LADA_LOADNSAVE_INITIALIZER_PARSE_VALUE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/utility/result_of.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_signed.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/not.hpp>
#include <boost/static_assert.hpp>

#include <boost/xpressive/match_results.hpp>
#include <boost/xpressive/regex_primitives.hpp>
#include <boost/xpressive/regex_algorithms.hpp>
#include <boost/xpressive/regex_actions.hpp>
#include <boost/xpressive/regex_compiler.hpp>

#include "../tree/tree.h"

#include "../grammar/grammar.h"
#include "convert_value_to_regex.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace initializer
    {

#     ifdef LADA_BODY
#       error LADA_BODY macro already exists.
#     endif
#     ifdef LADA_SIG
#       error LADA_SIG macro already exists.
#     endif
#     ifdef LADA_EXP
#       error LADA_EXP macro already exists.
#     endif
#     define LADA_SIG  bx::sregex const sig = (bx::set= '+', '-')
#     define LADA_EXP \
          bx::sregex const exp =    bx::as_xpr('.') >> !(bx::set='d','D','e','E')\
                                 >> !sig >> +bx::_d
#     define LADA_BODY( preregs, regex )\
        {                                  \
          namespace bx = boost::xpressive; \
          std::string value;               \
          preregs;                         \
          bx::sregex re = regex[bx::ref(value)=bx::s1]; \
          if( not bx::regex_match( _string, re ) ) return false; \
          _value = boost::lexical_cast<T_TYPE>( value ); \
          return true; \
        }

#     define BODY( regex )\
        
      //! Parses a signed integer.
      template<class T_TYPE> 
        typename boost::enable_if
        <
          typename boost::mpl::and_
          <
            typename boost::is_integral<T_TYPE> :: type,
            typename boost::is_signed<T_TYPE> :: type 
          > :: type, bool
        > :: type parse_value( std::string const& _string, T_TYPE &_value )
          { LADA_BODY( LADA_SIG;, (bx::s1 = (!sig >> +bx::_d)) ); }

      //! Parses an unsigned integer.
      template<class T_TYPE> 
        typename boost::enable_if
        <
          typename boost::mpl::and_
          <
            typename boost::is_integral<T_TYPE> :: type,
            typename boost::mpl::not_< typename boost::is_signed<T_TYPE> :: type > :: type
          > :: type, bool
        > :: type parse_value( std::string const& _string, T_TYPE &_value  )
         { LADA_BODY(;, (bx::s1 = (+bx::_d)) ); }

      //! Parses a floating point.
      template<class T_TYPE> 
        typename boost::enable_if
        <
          typename boost::is_floating_point<T_TYPE> :: type, bool
        > :: type parse_value( std::string const& _string, T_TYPE &_value  )
          { LADA_BODY( LADA_SIG; LADA_EXP;, (bx::s1 = (!sig >> +bx::_d >> !exp)) ); }

#     undef LADA_BODY
#     undef LADA_EXP
#     undef LADA_SIG


      //! Parses a string.
      inline bool parse_value( std::string const& _string, char const* _value  )
      {
        _value = _string.c_str();
        return true;
      }

      //! Parses a string.
      inline bool parse_value( std::string const& _string, std::string& _value  )
      {
        _value = _string;
        return true;
      }

      //! Converts a string to a bool
      inline bool parse_value( std::string const &_value, bool &_action)
      {
        namespace bx = boost::xpressive;
        bx::sregex const re = bx::sregex::compile( convert_value_to_regex( true ) );
        _action = bx::regex_match( _value, re );
        return true;
      }

      //! Metafunction to distinguish between actions and non-actions.
      template< class T, class ENABLE = void > struct is_action 
      {
        typedef boost::mpl::false_::type type;
      };
      //! Metafunction to distinguish between actions and non-actions.
      template< class T> struct is_action<T, typename T::action> 
      {
        typedef boost::mpl::true_::type type;
      };


      template< class T_TYPE >
        typename boost::enable_if
        <
          typename is_action<T_TYPE>::type, bool
        >::type parse_value( std::string const &_value, T_TYPE& _action )
          { return _action(_value); }
      
      //! Does nothing.
      inline bool parse_value( std::string const &_value,
                               boost::mpl::bool_<false> const &) { return true; }
    } // namespace initializer.
  } // namespace load_n_save

} // namespace LaDa


#endif
