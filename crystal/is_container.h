#ifndef LADA_CRYSTAL_IS_CONTAINER_H
#define LADA_CRYSTAL_IS_CONTAINER_H

#include "LaDaConfig.h"

#include <sstream>
#include <string>
#include <set>

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/not.hpp>
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

      //! Inherits from boost::mpl::true_ if template argument is a set.
      template<class T_CONTAINER>
        struct is_set : public boost::mpl::false_ {};
      //! Inherits from boost::mpl::true_ if template argument is a set.
      template<class K, class C, class A>
        struct is_set<std::set<K, C, A> >: public boost::mpl::true_ {};

      //! Inherits from boost::mpl::true_ if template argument is string convertible.
      template<class T> struct is_string : public boost::is_convertible<T, std::string> {};

      //! Inherits from boost::mpl::true_ if template argument is iterable
      template<class T> struct is_iterable 
        : public boost::mpl::or_< is_container<T>, is_set<T> > {};

      //! Inherits from boost::mpl::true_ if template argument is non-string scalar.
      template<class T> struct is_nonstring_scalar 
        : public boost::mpl::and_
                 < 
                   boost::mpl::not_< is_iterable<T> >, 
                   boost::mpl::not_< is_string<T> >
                 > {};

      //! Inherits from boost::mpl::true_ if template argument is scalar.
      template<class T> struct is_scalar 
        : public boost::mpl::not_< is_iterable<T> > {};

      //! \brief Helper function to print occupations.
      //! \details Specialization for scalar strings.
      template <class T>
        typename boost::enable_if<is_string<T>, std::string>::type
          print_occupation(T const &_in) { return _in; }
      //! \brief Helper function to print occupations.
      //! \details Specialization for scalar non-strings.
      template <class T>
        typename boost::enable_if<is_nonstring_scalar<T>, std::string>::type
          print_occupation(T const &_in) { std::ostringstream sstr; sstr << _in; return sstr.str(); }
      //! \brief Helper function to print occupations.
      //! \details Specialization for containers.
      template <class T>
        typename boost::enable_if<is_iterable<T>, std::string>::type
          print_occupation(T const &_in)
          {
            size_t const N(_in.size());
            if(N == 0) return "";
            std::ostringstream sstr;
            typename T::const_iterator i_first = _in.begin();
            typename T::const_iterator const i_end = _in.end();
            if(N == 1) return print_occupation(*i_first);
            sstr << *i_first;
            for(++i_first; i_first != i_end; ++i_first) sstr << ", " << *i_first;
            return sstr.str();
          }
    }

  } // namespace Crystal
} // namespace LaDa
  
#endif


