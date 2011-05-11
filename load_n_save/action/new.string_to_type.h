//
//  Version: $Id: string_to_type.h 1250 2009-07-26 21:04:07Z davezac $
//

#ifndef _LADA_LOADNSAVE_XPR_STRING_TO_TYPE_H_
#define _LADA_LOADNSAVE_XPR_STRING_TO_TYPE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/utility/enable_if.hpp>
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

#include "../string_type.h"

#include "type_to_regex.h"

namespace LaDa 
{
  namespace load_n_save
  {
    template< class T_TYPE, class T_BOOL = void > class StringToType;
#   ifdef LADA_BODY
#     error LADA_BODY macro already exists.
#   endif
#   ifdef LADA_SIG
#     error LADA_SIG macro already exists.
#   endif
#   ifdef LADA_EXP
#     error LADA_EXP macro already exists.
#   endif
#   define LADA_SIG  bx::sregex const sig = (bx::set= '+', '-')
#   define LADA_EXP \
        bx::sregex const exp =    bx::as_xpr('.') \
                               >> !( (bx::set='d','D','e','E') >> !sig)\
                               >> +bx::_d
#   define LADA_BODY( preregs, regex )\
      {                                  \
        namespace bx = boost::xpressive; \
        t_String value;                  \
        preregs;                         \
        bx::sregex re = regex[bx::ref(value)=bx::s1]; \
        if( not bx::regex_match( _string, re ) ) return false; \
        bx::sregex rep = (bx::set='d','D','E');  \
        std::string const reg( bx::regex_replace( _string, rep, std::string("e") ) ); \
        _value = boost::lexical_cast<T_TYPE>( reg ); \
        return true; \
      }

    //! Parses a signed integer.
    template<class T_TYPE> 
      class StringToType
            <
              T_TYPE,
              typename boost::enable_if
              <
                typename boost::mpl::and_
                <
                  typename boost::is_integral<T_TYPE> :: type,
                  typename boost::is_signed<T_TYPE> :: type 
                > :: type
              > :: type
            >
      {
        //! Functor.
        static bool apply( t_String const& _string, T_TYPE &_value )
          { LADA_BODY( LADA_SIG;, (bx::s1 = (!sig >> +bx::_d)) ); }
      };

    //! Parses an unsigned integer.
    template<class T_TYPE> 
      class StringToType
            <
              T_TYPE,
              typename boost::enable_if
              <
                typename boost::mpl::and_
                <
                  typename boost::is_integral<T_TYPE> :: type,
                  typename boost::mpl::not_< typename boost::is_signed<T_TYPE> :: type > :: type
                > :: type
              > :: type
            >
      {
        //! Functor.
        static bool apply( t_String const& _string, T_TYPE &_value )
          { LADA_BODY(;, (bx::s1 = (+bx::_d)) ); }
      };

    //! Parses a floating point.
    template<class T_TYPE> 
      class StringToType
            <
              T_TYPE,
              typename boost::enable_if
              <
                typename boost::is_floating_point<T_TYPE> :: type
              > :: type
            >
      {
        //! Functor.
        static bool apply( t_String const& _string, T_TYPE &_value )
          { LADA_BODY( LADA_SIG; LADA_EXP;, (bx::s1 = (!sig >> +bx::_d >> !exp)) ); }
      };

#   undef LADA_BODY
#   undef LADA_EXP
#   undef LADA_SIG


    //! Parses a string.
    template<> class StringToType<char const*, void>
    {
      //! Functor.
      static bool apply( t_String const& _string, char const* _value  )
      {
        _value = _string.c_str();
        return true;
      }
    };

    //! Parses a string.
    template<> class StringToType<t_String, void>
    {
      //! Functor.
      static bool apply( t_String const& _string, t_String &_value  )
      {
        _value = _string.c_str();
        return true;
      }
    };

    //! Converts a string to a bool
    template<> class StringToType<bool, void>
    {
      //! Functor.
      static bool apply( t_String const& _string, bool &_value  )
      {
        if( _string.size() == 0 ) { _value = true; return true; }
        namespace bx = boost::xpressive;
        bx::sregex const true_ = bx::sregex::compile( type_to_regex( true ) );
        bx::sregex const false_ = bx::sregex::compile( type_to_regex( false ) );
        if(bx::regex_match( _string, true_ ) ) _value = true;
        else if(bx::regex_match( _string, false_ ) ) _value = false;
        else return false;
        return true;
      }
    };

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
      class StringToType
            <
              T_TYPE,
              typename boost::enable_if
              <
                typename is_action<T_TYPE>::type
              > :: type
            >
      {
        //! Functor.
        static bool apply( t_String const& _string, T_TYPE &_action )
          { return _action(_string); }
      };
    
    //! Does nothing.
    template<> class StringToType<boost::mpl::bool_<false>, void>
    {
      //! Functor.
      static bool apply( t_String const&, boost::mpl::bool_<false> const &  )
        { return true; }
    };
  } // namespace load_n_save

} // namespace LaDa


#endif
