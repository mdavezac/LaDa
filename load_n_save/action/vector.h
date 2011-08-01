#ifndef LADA_LNS_VECTOR_H
#define LADA_LNS_VECTOR_H

#include <vector>

#include<boost/algorithm/string/split.hpp>

#include "string_to_type.h"
#include "type_to_string.h"
#include "type_to_regex.h"
#include "../string_type.h"

namespace LaDa
{
  namespace load_n_save
  {
    //! \brief Regex to parse a vector in an action*
    //! \details It is assumed that different elements are separated by spaces.
    //!          In other words, the elements themselves should not contain
    //!          spaces.
    template<class T_TYPE, class T_ALLOC>
      struct TypeToRegex<std::vector<T_TYPE, T_ALLOC>, void>
      {
        //! Returns regex string.
        static t_String apply()
          { return "(\\s*" + TypeToRegex<T_TYPE>::apply() + "\\s*)*"; }
      };

    //! \brief Parses a string into a type.
    //! \details It is assumed that different elements are separated by spaces.
    //!          In other words, the elements themselves should not contain
    //!          spaces.
    template<class T_TYPE, class T_ALLOC>
      struct StringToType<std::vector<T_TYPE, T_ALLOC>, void>
      {
        //! Functor.
        static bool apply( t_String const& _string, std::vector<T_TYPE, T_ALLOC> &_value )
        {
          namespace ba = boost::algorithm;
          std::vector<std::string> splitted;
          ba::split(splitted, _string, ba::is_any_of(" "), ba::token_compress_on);
          std::vector<std::string>::const_iterator i_first = splitted.begin();
          std::vector<std::string>::const_iterator const i_last = splitted.end();
          _value.resize(splitted.size());
          typename std::vector<T_TYPE, T_ALLOC>::iterator i_val = _value.begin();
          for(; i_first != i_last; ++i_first, ++i_val)
            StringToType<T_TYPE>::apply(*i_first, *i_val);
        }
      };

    //! \brief Dumps a vector to a string.
    //! \details Different elements are separated by spaces.
    //!          In other words, the elements themselves should not contain
    //!          spaces.
    template<class T_TYPE, class T_ALLOC>
      struct TypeToString<std::vector<T_TYPE, T_ALLOC>, void>
      {
        static t_String apply(std::vector<T_TYPE, T_ALLOC> const &_value)
        {
          if(_value.size() == 0) return "";
          typename std::vector<T_TYPE, T_ALLOC>::const_iterator i_first = _value.begin();
          typename std::vector<T_TYPE, T_ALLOC>::const_iterator const i_end = _value.end();
          t_String result = TypeToString<T_TYPE>::apply(*i_first);
          for(++i_first; i_first != i_end; ++i_first) result += " " + TypeToString<T_TYPE>::apply(*i_first);
          return result;
        }
      };
  } // namespace load_n_save
} // namespace LaDa

#endif 
