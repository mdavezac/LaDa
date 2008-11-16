//
//  Version: $Id$
//
#ifndef _LADA_DARWIN_OPERATORS_MAKE_POPULATOR_FUNCTOR_H_
#define _LADA_DARWIN_OPERATORS_MAKE_POPULATOR_FUNCTOR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/remove_const.hpp>

#include "discriminate.h"
#include "callback.h"

//! \cond
template< class EOT > class eoPopulator;
//! \endcond

namespace LaDa
{
  namespace GA
  {
    namespace Operator
    {
      //! Creates a genetic operator functor.
      template< class T_INDIVIDUAL, class T_POPULATOR = eoPopulator< T_INDIVIDUAL > >
        struct MakePopulator
        {
          //! A tag.
          template< size_t > class Tag {};
          public:
            //! Type of the individual.
            typedef T_INDIVIDUAL t_Individual;
            //! Type of the populator iterator.
            typedef T_POPULATOR t_Populator;
            //! The result function.
            typedef boost::function<void(t_Populator&)> t_PopFunctor;
           
            template< class T_FUNCTOR >
              static void transform( const T_FUNCTOR &_functor, t_PopFunctor &_result );
            //! Static functions which calls a populator functor with a populator.
            template< class T_FUNCTOR >
              static void call( const T_FUNCTOR &_functor, t_Populator &_populator );

            //! Unoverloaded function to call unary with populator.
            template< class T_FUNCTOR >
              void static call_unary( const T_FUNCTOR& _functor, T_POPULATOR &_populator );
            //! Unoverloaded function to call binary with populator.
            template< class T_FUNCTOR >
              void static call_binary( const T_FUNCTOR& _functor, T_POPULATOR &_populator );

            //! Unoverloaded function to transform unary functor.
            template< class T_FUNCTOR >
              void static transform_unary( const T_FUNCTOR &_functor, t_PopFunctor &_result );
            //! Unoverloaded function to call binary with populator.
            template< class T_FUNCTOR >
              void static transform_binary( const T_FUNCTOR& _functor, t_PopFunctor &_result );
          private:
            //! A unary tag type.
            typedef Tag< Disc::unary_index > t_unary;
            //! A binary tag type.
            typedef Tag< Disc::binary_index > t_binary;
            //! A const unary tag type.
            typedef Tag< Disc::const_unary_index > t_const_unary;
            //! A const binary tag type.
            typedef Tag< Disc::const_binary_index > t_const_binary;
            //! An "other" tag type.
            typedef Tag< Disc::other_index > t_other;
            //! Static functions which transforms a unary functor into a populator.
            template< class T_FUNCTOR >
              static void transform_( const T_FUNCTOR &_functor,
                                      t_PopFunctor &_result, 
                                      const t_unary& );
            //! Static functions which transforms a binary functor into a populator.
            template< class T_FUNCTOR >
              static void transform_( const T_FUNCTOR &_functor,
                                      t_PopFunctor &_result, 
                                      const t_binary& );
            //! Static functions which transforms a const unary functor into a populator.
            template< class T_FUNCTOR >
              static void transform_( const T_FUNCTOR &_functor,
                                      t_PopFunctor &_result, 
                                      const t_const_unary& );
            //! Static functions which transforms a const binary functor into a populator.
            template< class T_FUNCTOR >
              static void transform_( const T_FUNCTOR &_functor,
                                      t_PopFunctor &_result, 
                                      const t_const_binary& );
            //! Static functions which transforms a populator functor into a populator functor.
            template< class T_FUNCTOR >
              static void transform_( const T_FUNCTOR &_functor,
                                      t_PopFunctor &_result, 
                                      const t_other& );

            //! Calls a unary functor with a populator
            template< class T_FUNCTOR >
              static void call_( const T_FUNCTOR& _functor, T_POPULATOR &_populator, const t_unary& );
            //! Calls a binary functor with a populator
            template< class T_FUNCTOR >
              static void call_( const T_FUNCTOR& _functor, T_POPULATOR &_populator, const t_binary& );
            //! Calls a const unary functor with a populator
            template< class T_FUNCTOR >
              static void call_( const T_FUNCTOR& _functor,
                                 T_POPULATOR &_populator, 
                                 const t_const_unary& );
            //! Calls a const binary functor with a populator
            template< class T_FUNCTOR >
              static void call_( const T_FUNCTOR& _functor,
                                 T_POPULATOR &_populator, 
                                 const t_const_binary& );
            //! Calls a populator functor with a populator
            template< class T_FUNCTOR >
              static void call_( const T_FUNCTOR& _functor, T_POPULATOR &_populator, const t_other& );
        };

#     ifdef INPOP
#       error INPOP already defined.
#     endif
#     define INPOP \
      template< class T_INDIVIDUAL, class T_POPULATOR > template< class T_FUNCTOR > \
        void MakePopulator<T_INDIVIDUAL, T_POPULATOR>  

      INPOP :: transform_unary( const T_FUNCTOR &_functor, t_PopFunctor &_result )
        { _result = constUnaryCallBack<t_Individual, t_Populator>( _functor ); }
      INPOP :: transform_binary( const T_FUNCTOR &_functor, t_PopFunctor &_result )
        { _result = constBinaryCallBack<t_Individual, t_Populator>( _functor ); }
      INPOP :: call_unary( const T_FUNCTOR &_functor, t_Populator& _pop )
        { (constUnaryCallBack<t_Individual, t_Populator>( _functor ))( _pop ); }
      INPOP :: call_binary( const T_FUNCTOR &_functor, t_Populator& _pop  )
        { (constBinaryCallBack<t_Individual, t_Populator>( _functor ))( _pop ); }

      INPOP :: transform( const T_FUNCTOR &_functor, t_PopFunctor &_result )
        { transform_( _functor, _result, Tag<Discriminate<T_INDIVIDUAL, T_FUNCTOR >::value>() ); }
      INPOP :: call( const T_FUNCTOR &_functor, t_Populator &_populator )
        { call_( _functor, _populator, Tag<Discriminate<T_INDIVIDUAL, T_FUNCTOR >::value>() ); }

      INPOP :: transform_( const T_FUNCTOR &_functor, t_PopFunctor &_result, const t_unary& )
        { transform_unary( _functor, _result ); }
      INPOP :: transform_( const T_FUNCTOR &_functor, t_PopFunctor &_result, const t_binary& )
        { transform_binary( _functor, _result ); }
      INPOP :: transform_( const T_FUNCTOR &_functor, t_PopFunctor &_result, const t_const_unary& )
        { transform_unary( _functor, _result ); }
      INPOP :: transform_( const T_FUNCTOR &_functor,
                           t_PopFunctor &_result,
                           const t_const_binary& )
        { transform_binary( _functor, _result ); }
      INPOP :: transform_( const T_FUNCTOR &_functor, t_PopFunctor &_result, const t_other& )
        { _result = _functor; }

      INPOP :: call_( const T_FUNCTOR &_functor, t_Populator& _pop, const t_unary& )
        { call_unary( _functor, _pop ); }
      INPOP :: call_( const T_FUNCTOR &_functor, t_Populator& _pop, const t_binary& )
        { call_binary( _functor, _pop ); }
      INPOP :: call_( const T_FUNCTOR &_functor, t_Populator& _pop, const t_const_unary& )
        { call_unary( _functor, _pop ); }
      INPOP :: call_( const T_FUNCTOR &_functor, t_Populator& _pop, const t_const_binary& )
        { call_binary( _functor, _pop ); }
      INPOP :: call_( const T_FUNCTOR &_functor, t_Populator& _pop, const t_other& )
        { _functor( _pop ); }


#     undef INPOP
    } // namespace Operators

  } // namespace GA
} // namespace LaDa



#endif 
