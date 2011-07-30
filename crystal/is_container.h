#ifndef LADA_CRYSTAL_IS_CONTAINER_H
#define LADA_CRYSTAL_IS_CONTAINER_H

#include "LaDaConfig.h"

#include <sstream>

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/bool.hpp>

namespace LaDa
{
  namespace crystal
  {
    namespace details
    {
      //! Inherits from boost::mpl::true_ if template argument is a container.
      template<class T_CONTAINER>
        struct is_container : public boost::mpl::false_ {};
      //! Inherits from boost::mpl::true_ if template argument is a container.
      template< template<class A, class B> class T_CONTAINER, class TYPE, class ALLOC>
        struct is_container<T_CONTAINER<TYPE, ALLOC> >: public boost::mpl::true_ {};

      template <class T>
        typename boost::disable_if<is_container<T>, std::string>::type
          print_container(T const &_in) { return std::ostringstream(_in).str(); }
      template <class T>
        typename boost::enable_if<is_container<T>, std::string>::type
          print_container(T const &_in)
          {
            if(_in.size() == 0) return "";
            if(_in.size() == 1) return std::ostringstream(_in.front()).str();
            std::ostringstream sstr;
            sstr << _in.front();
            typename T::const_iterator i_first = _in.begin();
            typename T::const_iterator const i_end = _in.end();
            for(++i_first; i_first != i_end; ++i_first) sstr << ", " << *i_first;
            return sstr.str();
          }
    }

  } // namespace Crystal
} // namespace LaDa
  
#endif


