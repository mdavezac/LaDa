#ifndef LADA_CRYSTAL_IS_CONTAINER_H
#define LADA_CRYSTAL_IS_CONTAINER_H

#include "LaDaConfig.h"

#include <sstream>
#include <string>

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/or.hpp>
#include <boost/type_traits/is_convertible.hpp>

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

      //! \brief Helper function to print occupations.
      //! \details Specialization for scalar strings.
      template <class T>
        typename boost::enable_if
          < 
             boost::mpl::and_
             < 
               boost::mpl::not_< is_container<T> >, 
               boost::is_convertible<T, std::string> 
             >, std::string
          > :: type
          print_occupation(T const &_in) { return _in; }
      //! \brief Helper function to print occupations.
      //! \details Specialization for scalar non-strings.
      template <class T>
        typename boost::enable_if
        < 
          boost::mpl::and_
          < 
            boost::mpl::not_< is_container<T> >,
            boost::mpl::not_<boost::is_convertible<T, std::string> > 
          >, std::string
        > :: type
          print_occupation(T const &_in) { std::ostringstream sstr; sstr << _in; return sstr.str(); }
      //! \brief Helper function to print occupations.
      //! \details Specialization for containers.
      template <class T>
        typename boost::enable_if<is_container<T>, std::string>::type
          print_occupation(T const &_in)
          {
            if(_in.size() == 0) return "";
            if(_in.size() == 1) return print_occupation(_in.front());
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


