#ifndef LADA_LOADNSAVE_BITWISE_MAP_ACTION_H
#define LADA_LOADNSAVE_BITWISE_MAP_ACTION_H

#include "LaDaConfig.h"

#include <cmath>
#include <map>
#include <algorithm>
#include <iterator>
#include <utility>

#include <boost/xpressive/regex_compiler.hpp>
#include <boost/xpressive/regex_algorithms.hpp>

#include "../string_type.h"
#include "action_base.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace action_
    {
      //! \brief Enum-like action.
      //! \details The option can take the value of any key or value in the map.
      //!          If inclusive, the result is the compound of all values found
      //!          in the option, linked via |. If not inclusive, the option
      //!          should be a single value or key.
      template< class T_TYPE >
        class Enum
        {
          public:
            //! A trait wich says this is an action.
            typedef void action;
            //! The map type
            typedef std::map<t_String, T_TYPE> t_map;
            //! \brief Constructor. 
            //! \params _var
            Enum   ( T_TYPE &_var, t_map const& _map, bool _incl = true )
                 : var_(_var), map_(new t_map(_map)), inclusive_(_incl) {} 
            //! CopyConstructor. 
            Enum( Enum const &_c) : var_(_c.var_), map_(_c.map_), inclusive_(_c.inclusive_) {} 
            //! Parses string into variable.
            bool operator()( t_String const &_x ) const;
            //! Returns a regex string.
            t_String operator()() const;
            //! Assigns the default value.
            bool assign_default( t_String const& _default ) const
              { return operator()( _default ); }
            //! Assigns the default value.
            template<class T>
              bool assign_default( T const& _default ) const
                { var_ = _default; return true; }

            //! Prints a value
            virtual t_String str() const
              { return str(typename boost::is_unsigned<T_TYPE>::type()); }
     
          protected:
            //! \brief specialization for unsigned integers.
            //! \details Makes sure that in the case of inclusive variables, we
            //!          start from the smallest and end up with the largest
            //!          integer.  This way, compound statements can be
            //!          correctly detected, X|Y == "xy", whithout also
            //!          including compound statements such as "all" if "all" == "xyz". 
            t_String str(boost::mpl::bool_<true>) const;
            //! general case
            t_String str(boost::mpl::bool_<false>) const;
            //! Holds reference to variable.
            T_TYPE &var_;
            //! Holds a constant reference to the map.
            boost::shared_ptr<t_map> map_;
            //! Whether should match one or many items of the map
            bool inclusive_;
        };
     
      template<class T_TYPE>
        bool Enum<T_TYPE>::operator()( t_String const &_x ) const
        { 
          var_ = T_TYPE(0);
          namespace bx = boost::xpressive;
          typename t_map::const_iterator i_first = map_->begin();
          typename t_map::const_iterator const i_end = map_->end();
          bx::sregex re;
          bool found(false);
          for(; i_first != i_end and (inclusive_ or not found); ++i_first)
          {
            re = bx::sregex::compile(i_first->first);
            if( not bx::regex_search( _x, re ) ) continue;
            
            var_ |= i_first->second;
            found = true;
          }
          return found; 
        }
     
      template<class T_TYPE>
        t_String Enum<T_TYPE>::operator()() const
        {
          t_String result;
          typename t_map::const_iterator i_first = map_->begin();
          typename t_map::const_iterator const i_end = map_->end();
          if( i_first == i_end ) return "";
          result +=   "(" + i_first->first;
          for(++i_first; i_first != i_end; ++i_first)
            result +=   "|" + i_first->first;
          result +=   ")";
          if( inclusive_ ) result += "+";
          return result;
        }

      template<class T_TYPE>
        t_String Enum<T_TYPE>::str(boost::mpl::bool_<false>) const
        {
          typename t_map::const_iterator i_first = map_->begin();
          typename t_map::const_iterator const i_end = map_->end();
          if( i_first == i_end ) return "";
          // first checks that one variable does not fit the value.
          for(; i_first != i_end; ++i_first)
            if(var_ == i_first->second) return i_first->first;
          if(not inclusive_)
            BOOST_THROW_EXCEPTION( error::enum_transcript_error() 
                                       << error::option_name("Cannot parse enum input.") );
          t_String result;
          for(i_first = map_->begin(); i_first != i_end; ++i_first)
            if(var_ & i_first->second) result += i_first->first;
          return result;
        }

      namespace details
      {
        template<class T>
          bool sort_second( typename T::value_type const &_a,
                            typename T::value_type const &_b  ) 
            {return _a.second < _b.second; }
      };
      template<class T_TYPE>
        t_String Enum<T_TYPE>::str(boost::mpl::bool_<true>) const
        {
          typename t_map::const_iterator i_first = map_->begin();
          typename t_map::const_iterator const i_end = map_->end();
          T_TYPE var(var_);
          if( i_first == i_end ) return "";

          // first checks that one variable does not fit the value.
          for(; i_first != i_end; ++i_first)
            if(var == i_first->second) return i_first->first;
          if(not inclusive_)
            BOOST_THROW_EXCEPTION( error::enum_transcript_error() 
                                       << error::option_name("Cannot parse enum input.") );

          // Sort key-values according to absolute value of the value.
          typedef typename std::pair<typename t_map::key_type, typename t_map::mapped_type> vt;
          std::vector<vt> sorted(map_->size());
          typename std::vector<vt>::iterator i_sorted = sorted.begin();
          typename std::vector<vt>::iterator i_sorted_end = sorted.begin();
          for(i_first = map_->begin(); i_first != i_end; ++i_first, ++i_sorted)
          {
            i_sorted->first = i_first->first;
            i_sorted->second = i_first->second;
          }
          std::sort(sorted.begin(), sorted.end(), &details::sort_second<t_map>); 
          t_String result("");
          i_sorted = sorted.begin();
          i_sorted_end = sorted.end();
          for(; i_sorted != i_sorted_end; ++i_sorted)
            if(var & i_sorted->second)
            {
              var -= i_sorted->second;
              result += i_sorted->first;
            }
          return result;
        }
    } // namespace action.
    
    //! Returns an Enum action.
    template< class T_TYPE >
      action_::Enum<T_TYPE> enum_( T_TYPE &_a, 
                                   typename action_::Enum<T_TYPE>::t_map const& _map, 
                                   bool _inc = true )
        { return action_::Enum<T_TYPE>(_a, _map, _inc); }

  }
}

#endif
