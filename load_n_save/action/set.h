#ifndef LADA_LNS_ACTION_SET_H
#define LADA_LNS_ACTION_SET_H

#include <vector>
#include <set>

#include<boost/algorithm/string/classification.hpp>
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
    template<class KEY, class COMPARE, class ALLOC>
      struct TypeToRegex<std::set<KEY, COMPARE, ALLOC>, void>
      {
        //! Returns regex string.
        static t_String apply()
          { return "(\\s*" + TypeToRegex<KEY>::apply() + "\\s*)*"; }
      };

    //! \brief Parses a string into a type.
    //! \details It is assumed that different elements are separated by spaces.
    //!          In other words, the elements themselves should not contain
    //!          spaces.
    template<class KEY, class COMPARE, class ALLOC>
      struct StringToType<std::set<KEY, COMPARE, ALLOC>, void>
      {
        //! Functor.
        static bool apply( t_String const& _string, std::set<KEY, COMPARE, ALLOC> &_value )
        {
          namespace ba = boost::algorithm;
          std::vector<std::string> splitted;
          ba::split(splitted, _string, ba::is_any_of(" "), ba::token_compress_on);
          std::vector<std::string>::const_iterator i_first = splitted.begin();
          std::vector<std::string>::const_iterator const i_last = splitted.end();
          _value.clear();
          for(; i_first != i_last; ++i_first)
          {
            KEY key;
            StringToType<KEY>::apply(*i_first, key);
            _value.insert(key);
          }
        }
      };

    //! \brief Dumps a vector to a string.
    //! \details Different elements are separated by spaces.
    //!          In other words, the elements themselves should not contain
    //!          spaces.
    template<class KEY, class COMPARE, class ALLOC>
      struct TypeToString<std::set<KEY, COMPARE, ALLOC>, void>
      {
        static t_String apply(std::set<KEY, COMPARE, ALLOC> const &_value)
        {
          if(_value.size() == 0) return "";
          typename std::set<KEY, COMPARE, ALLOC>::const_iterator i_first = _value.begin();
          typename std::set<KEY, COMPARE, ALLOC>::const_iterator const i_end = _value.end();
          t_String result = TypeToString<KEY>::apply(*i_first);
          for(++i_first; i_first != i_end; ++i_first)
            result += " " + TypeToString<KEY>::apply(*i_first);
          return result;
        }
      };
  } // namespace load_n_save
} // namespace LaDa

#endif 
