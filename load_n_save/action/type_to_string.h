#ifndef LADA_LOADNSAVE_XPR_TYPE_TO_STRING_H
#define LADA_LOADNSAVE_XPR_TYPE_TO_STRING_H

#include "LaDaConfig.h"

#include "../string_type.h"

#include "type_to_regex.h"

namespace LaDa 
{
  namespace load_n_save
  {
    template< class T_TYPE, class T_BOOL = void > struct TypeToString;

    //! Parses a signed integer.
    template<class T_TYPE> 
      struct TypeToString<T_TYPE, void>
      {
        static t_String apply(T_TYPE const &_value)
        {
          std::ostringstream sstr;
          sstr << _value;
          return sstr.str();
        }
      };
    template<> struct TypeToString<t_String, void> 
    {
      static t_String apply(t_String const &_value) { return _value; }
    };
  } // namespace load_n_save

} // namespace LaDa


#endif
