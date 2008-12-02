//
//  Version: $Id$
//
#ifndef  _LADA_GA_STATIC_TRANSLATE_H_
#define  _LADA_GA_STATIC_TRANSLATE_H_

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <string>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/utility/enable_if.hpp>

namespace LaDa
{
  namespace Crystal
  {
    class Structure;
  }

  namespace GA
  {
    //! \cond
    namespace details
    {
      class isStaticTranslateDerived {};
    }
    //! \endcond


    //! \brief Adds operator<< to an object using a static translate.
    //! \details Objects created with this metafunction will call the correct
    //!          static translate class.
    template< class T_OBJECT, template< class > class T_STATICTRANSLATE >
      class static_translate : public T_OBJECT
      {
        //! Holds implementation details.
        struct implementation : public T_OBJECT, public details::isStaticTranslateDerived 
        {
          //! Type of the Translate.
          typedef T_STATICTRANSLATE< T_OBJECT > t_Translate;
        };

        public:
         //! New object type.
         typedef implementation type;
      };

    //! Overloads from object to string.
    template< class T_OBJECT >
      typename boost::enable_if
      <
        boost::is_base_of< details::isStaticTranslateDerived, T_OBJECT >
      > :: type operator<<( std::string &_str, const T_OBJECT& _o )
        { T_OBJECT :: t_Translate :: translate( _o, _str ); }
    //! Overloads from string to object.
    template< class T_OBJECT >
      typename boost::enable_if
      <
        boost::is_base_of< details::isStaticTranslateDerived, T_OBJECT >
      > :: type operator<<( T_OBJECT& _o,  const std::string &_str )
        { T_OBJECT :: t_Translate :: translate( _str, _o ); }
    //! Overloads from structure to object.
    template< class T_OBJECT >
      typename boost::enable_if
      <
        boost::is_base_of< details::isStaticTranslateDerived, T_OBJECT >
      > :: type operator<<( T_OBJECT& _o,  const Crystal::Structure &_str )
        { T_OBJECT :: t_Translate :: translate( _str, _o ); }
    //! Overloads from object to structure.
    template< class T_OBJECT >
      typename boost::enable_if
      <
        boost::is_base_of< details::isStaticTranslateDerived, T_OBJECT >
      > :: type operator<<( Crystal::Structure &_str, const T_OBJECT& _o )
        { T_OBJECT :: t_Translate :: translate( _o, _str ); }
  }
} // namespace LaDa
#endif

