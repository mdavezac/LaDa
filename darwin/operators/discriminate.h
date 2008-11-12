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
           //! Can be instantiated by Unary operators.
           template<class U, bool (T_FUNCTOR::*)( T_INDIVIDUAL& ) const> struct cUnary {};
           //! Can be instantiated by Binary operators.
           template<class U, bool (T_FUNCTOR::*)( T_INDIVIDUAL&, const T_INDIVIDUAL& ) const>
             struct cBinary {};

           //! \cond
           typedef char t_other;
           union t_unary { char i[2]; };
           union t_binary { char i[3]; };
           union t_const_unary { char i[4]; };
           union t_const_binary { char i[5]; };
           //! \endcond

           //! This overload can be instantiated with any class.
           template<class U> static t_other Test(...);
           //! This overload can be instantiated only with the correct class and class member.
           template<class U> static t_unary Test( Unary<U, &U::operator() >*);
           //! This overload can be instantiated only with the correct class and class member.
           template<class U> static t_binary Test( Binary<U, &U::operator() >*);
           //! This overload can be instantiated only with the correct class and class member.
           template<class U> static t_const_unary Test( cUnary<U, &U::operator() >*);
           //! This overload can be instantiated only with the correct class and class member.
           template<class U> static t_const_binary Test( cBinary<U, &U::operator() >*);
           public:
             //! The resulting value.
             static const size_t value = sizeof(Test<T_FUNCTOR>(0));
             //! True if the functor is other
             static const bool is_other = sizeof(Test<T_FUNCTOR>(0)) == sizeof( t_other );
             //! True if the functor is unary
             static const bool is_unary = sizeof(Test<T_FUNCTOR>(0)) == sizeof( t_unary );
             //! True if the functor is binary
             static const bool is_binary = sizeof(Test<T_FUNCTOR>(0)) == sizeof( t_binary );
             //! True if the functor is const unary
             static const bool is_const_unary =    sizeof(Test<T_FUNCTOR>(0))
                                                == sizeof( t_const_unary );
             //! True if the functor is const binary
             static const bool is_const_binary =    sizeof(Test<T_FUNCTOR>(0))
                                                 == sizeof( t_const_binary );
             //! Index of other functors.
             static const size_t other_index = sizeof( t_other );
             //! Index of unary functors.
             static const size_t unary_index = sizeof( t_unary );
             //! Index of binary functors.
             static const size_t binary_index = sizeof( t_binary );
             //! Index of const unary functors.
             static const size_t const_unary_index = sizeof( t_const_unary );
             //! Index of const binary functors.
             static const size_t const_binary_index = sizeof( t_const_binary );
        };
     
      //! A conditional metafunction which is true if T_FUNCTOR is Unary.
      template< class T_INDIVIDUAL, class T_FUNCTOR >
        struct IsUnary
        {
          //! result.
          const static bool value = Discriminate<T_INDIVIDUAL, T_FUNCTOR>::is_unary;
        };
      //! A conditional metafunction which is true if T_FUNCTOR is Binary.
      template< class T_INDIVIDUAL, class T_FUNCTOR >
        struct IsBinary
        {
          //! result.
          const static bool value = Discriminate<T_INDIVIDUAL, T_FUNCTOR>::is_binary;
        };
      //! A conditional metafunction which is true if T_FUNCTOR is Binary.
      template< class T_INDIVIDUAL, class T_FUNCTOR >
        struct IsOther
        {
          //! result.
          const static bool value = Discriminate<T_INDIVIDUAL, T_FUNCTOR>::is_other;
        };
    }

  }
}



#endif 
