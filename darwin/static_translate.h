//
//  Version: $Id$
//
#ifndef  _LADA_GA_STATIC_TRANSLATE_H_
#define  _LADA_GA_STATIC_TRANSLATE_H_

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <string>

namespace LaDa
{
  namespace Crystal
  {
    class Structure;
  }

  namespace GA
  {
    //! \brief Adds operator<< to an object using a static translate.
    //! \details Objects created with this metafunction will call the correct
    //!          static translate class.
    template< class T_OBJECT, template< class > class T_STATICTRANSLATE >
      class static_translate : public T_OBJECT
      {
        //! Holds implementation details.
        class implementation : public T_OBJECT {};

        public:
         //! Type of the object.
         typedef T_OBJECT t_Object;
         //! Type of the Translate.
         typedef T_STATICTRANSLATE< t_Object > t_Translate;

         //! New object type.
         typedef implementation type;
      };

    //! Overloads from object to string.
    template< class T_OBJECT, template< class > class T_ST >
      void operator<<( std::string &_str,
                       const typename static_translate< T_OBJECT, T_ST > :: type& _o )
        { static_translate< T_OBJECT, T_ST > :: t_Translate :: translate( _o, _str ); }
    //! Overloads from string to object.
    template< class T_OBJECT, template< class > class T_ST >
      void operator<<( typename static_translate< T_OBJECT, T_ST > :: type& _o, 
                       std::string &_str )
        { static_translate< T_OBJECT, T_ST > :: t_Translate :: translate( _str, _o ); }
    //! Overloads from structure to object.
    template< class T_OBJECT, template< class > class T_ST >
      void operator<<( Crystal::Structure &_str,
                       const typename static_translate< T_OBJECT, T_ST > :: type& _o )
        { static_translate< T_OBJECT, T_ST > :: t_Translate :: translate( _str, _o ); }
    //! Overloads from object to structure.
    template< class T_OBJECT, template< class > class T_ST >
      void operator<<( typename static_translate< T_OBJECT, T_ST > :: type& _o, 
                       const Crystal::Structure& _str )
        { static_translate< T_OBJECT, T_ST > :: t_Translate :: translate( _o, _str ); }
  }
} // namespace LaDa
#endif

