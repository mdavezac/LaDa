//
//  Version: $Id: type_to_regex.h 1293 2009-09-08 05:51:37Z davezac $
//

#ifndef _LADA_LOADNSAVE_XPR_TYPE_TO_REGEX_H_
#define _LADA_LOADNSAVE_XPR_TYPE_TO_REGEX_H_


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
#include <boost/lexical_cast.hpp>
#include "../string_type.h"


namespace LaDa 
{
  namespace load_n_save
  {
    //! Transforms a type into a regex string.
    template<class T_TYPE, class T_VOID=void> struct TypeToRegex;
    //! Regex for signed integers.
    template<class T_TYPE>
      struct TypeToRegex
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
        //! Returns regex string.
        static t_String apply() { return  "[\\+-]?\\d+"; }
      };

    //! Regex for unsigned integers.
    template<class T_TYPE> 
      struct TypeToRegex
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
        //! Returns regex string.
        static t_String apply()   { return  "\\d+"; }
      };

    //! Regex for floating points.
    template<class T_TYPE> 
      struct TypeToRegex
      <
        T_TYPE,
        typename boost::enable_if
        <
          typename boost::is_floating_point<T_TYPE> :: type
        > :: type
      >
      {
        //! Returns regex string.
        static t_String apply() { return  "[\\+-]?\\d+(\\.([dDeE][\\+-]?)?\\d+)?"; }
      };
       
    
    //! Regex for a string.
    template<> struct TypeToRegex<t_String, void> 
    {
      //! Returns regex string.
      static t_String apply() { return ".*"; }
    };
      
    //! Regex for a string.
    template<> struct TypeToRegex<char const*, void> 
    {
      //! Returns regex string.
      static t_String apply() { return ".*"; }
    };

    //! Regex for a boolean.
    template<> struct TypeToRegex<bool, void> 
    {
      //! Returns regex string.
      static t_String apply() { return "((?i:(true|t))|[1-9]+|(?i:(false|f))|0+)?"; }
    };
    
  }
}

#endif
