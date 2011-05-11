//
//  Version: $Id: convert_value_to_regex.h 1189 2009-06-17 02:19:52Z davezac $
//

#ifndef _LADA_LOADNSAVE_INITIALIZER_CONVERT_VALUE_TO_REGEX_H_
#define _LADA_LOADNSAVE_INITIALIZER_CONVERT_VALUE_TO_REGEX_H_


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/lexical_cast.hpp>


namespace LaDa 
{
  namespace load_n_save
  {
    namespace initializer
    {
      //! Converts a value to a string via boost::lexical_cast.
      template< class T_VALUE >
        std::string convert_value_to_regex( T_VALUE const &_value ) 
          { return boost::lexical_cast< std::string >( _value ); }
      //! Returns \a _value.
      inline std::string convert_value_to_regex( std::string const &_value ) 
        { return _value; }
      //! Returns \a _value.
      inline std::string convert_value_to_regex( char const *_value ) 
        { return _value; }
      //! Converts a bool to a regex.
      inline std::string convert_value_to_regex( bool _value ) 
        { return _value ? "(?i:(true|t))|[1-9]+": "(?i:(false|f))|0+"; }
    }
  }
}

#endif
