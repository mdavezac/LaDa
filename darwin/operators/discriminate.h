//
//  Version: $Id$
//
#ifndef _LADA_DARWIN_OPERATORS_GENETIC_H_
#define _LADA_DARWIN_OPERATORS_GENETIC_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LaDa
{
  namespace GA
  {

    namespace Operator
    {
      //! Metafunction to distinguish between operator types.
      template< class T_INDIVIDUAL, class T_FUNCTOR >
        class Discriminate
        {
           //! Can be instantiated by Unary operators.
           template<class U, bool (T_FUNCTOR::*)( T_INDIVIDUAL& )> struct Unary {};
           //! Can be instantiated by Binary operators.
           template<class U, bool (T_FUNCTOR::*)( T_INDIVIDUAL&, const T_INDIVIDUAL& )>
             struct Binary {};
           //! This overload can be instantiated with any class.
           template<class U> static char Test(...);
           //! This overload can be instantiated only with the correct class and class member.
           template<class U> static char[2] Test( Unary<U, &U::operator() >*);
           //! This overload can be instantiated only with the correct class and class member.
           template<class U> static char[3] Test( Binary<U, &U::operator() >*);
           public:
             //! The resulting value.
             static const int value = sizeof(Test<T_FUNCTOR>(NULL))
        };
     
      //! A conditional metafunction which is true if T_FUNCTOR is Unary.
      template< class T_INDIVIDUAL, class T_FUNCTOR >
        struct IsUnary
        {
          //! result.
          const static bool value = FunctorType == 2;
        };
      //! A conditional metafunction which is true if T_FUNCTOR is Binary.
      template< class T_INDIVIDUAL, class T_FUNCTOR >
        struct IsBinary
        {
          //! result.
          const static bool value = FunctorType == 3;
        };
    }

  }
}



#endif 
