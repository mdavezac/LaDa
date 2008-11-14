//
//  Version: $Id$
//
#ifndef _LADA_DARWIN_OPERATORS_POPULATOR_FUNCTOR_H_
#define _LADA_DARWIN_OPERATORS_POPULATOR_FUNCTOR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/remove_const.hpp>

#include <opt/modifiers.h>

#include "discriminate.h"

//! \cond
template< class EOT > class eoPopulator;
//! \endcond


namespace LaDa
{
  namespace GA
  {
    //! Holds general stuff for GA operators.
    namespace Operator
    {
      //! \cond
      namespace details
      {
        template< size_t D > struct branch {};
      }
      //! \endcond

      template< class T_DERIVED, class T_INDIVIDUAL, class T_POPULATOR > class PopulatorFunctor;

      //! Assigns a callback functor to a boost function.
      template< class T_DERIVED, class T_INDIVIDUAL, class T_POPULATOR >
        void assign( PopulatorFunctor<T_DERIVED, T_INDIVIDUAL, T_POPULATOR>& _b, 
                     boost::function< void( T_POPULATOR& ) > &_a );
      //! Assigns a const callback functor to a boost function.
      template< class T_DERIVED, class T_INDIVIDUAL, class T_POPULATOR >
        void assign( const PopulatorFunctor<T_DERIVED, T_INDIVIDUAL, T_POPULATOR>& _b, 
                     boost::function< void( T_POPULATOR& ) > &_a );


      //! Creates a genetic operator, eg, a functor of a T_POPULATOR
      template< class T_DERIVED,
                class T_INDIVIDUAL = typename T_DERIVED :: t_Individual,
                class T_POPULATOR = eoPopulator< T_INDIVIDUAL > >
        class PopulatorFunctor
        {
            friend void assign<T_DERIVED, T_INDIVIDUAL, T_POPULATOR>
               ( PopulatorFunctor<T_DERIVED, T_INDIVIDUAL, T_POPULATOR>&,
                 boost::function<void(T_POPULATOR&)>& );
            friend void assign<T_DERIVED, T_INDIVIDUAL, T_POPULATOR>
               ( const PopulatorFunctor<T_DERIVED, T_INDIVIDUAL, T_POPULATOR>&,
                 boost::function<void(T_POPULATOR&)>& );
          private:
            //! Type of the discrimination metafunction.
            typedef Discriminate<T_INDIVIDUAL, T_DERIVED> t_Discriminate;
            //! "Discrimination" index of the functor.
            const static size_t discrimination = t_Discriminate :: value;
            //! A class to distinguish between types.
            template< size_t D > struct branch;
     
          public:
            //! Type of the derived class.
            typedef T_DERIVED t_Derived;
            //! Type of the indidividual
            typedef T_INDIVIDUAL t_Individual;
            //! Type of the Populator
            typedef T_POPULATOR t_Populator;
     
            //! Virtual destrcutor, just to make sure.
            virtual ~PopulatorFunctor() {}
     
            //! Functor over populator. Branches to correct format for t_Derived.
            void operator()( t_Populator& _populator )
              { details::branch< discrimination >::call( this, _populator ); }
            //! Functor over populator. Branches to correct format for t_Derived.
            void operator()( t_Populator& _populator ) const
              { details::branch< discrimination > :: call( this, _populator ); }
            
          private:
            //! Assigns this object to \a _function.
            void assign( boost::function< void( t_Populator& ) >& );
        };

      template< class T_DERIVED, class T_INDIVIDUAL, class T_POPULATOR >
        void assign( PopulatorFunctor<T_DERIVED, T_INDIVIDUAL, T_POPULATOR>& _b, 
                     boost::function< void( T_POPULATOR& ) > &_a )
        {
          typedef PopulatorFunctor< T_DERIVED, T_INDIVIDUAL, T_POPULATOR > t_PopFunc;
          const T_DERIVED * const derived = static_cast< T_DERIVED* >( &_b );
          void (t_PopFunc::*ptr_popop)( T_POPULATOR& ) = &t_PopFunc :: operator();
          _a = boost::bind( ptr_popop, *derived, _1 );
        }

      template< class T_DERIVED, class T_INDIVIDUAL, class T_POPULATOR >
        void assign( const PopulatorFunctor<T_DERIVED, T_INDIVIDUAL, T_POPULATOR>& _b, 
                     boost::function< void( T_POPULATOR& ) > &_a )
        {
          typedef PopulatorFunctor< T_DERIVED, T_INDIVIDUAL, T_POPULATOR > t_PopFunc;
          const T_DERIVED * const derived = static_cast< T_DERIVED* >( &_b );
          void (t_PopFunc::*ptr_popop)( T_POPULATOR& ) const = &t_PopFunc :: operator();
          _a = boost::bind( ptr_popop, *derived, _1 );
        }

      namespace details
      {
        template<> struct branch< Disc::unary_index > 
        {
          template< class T_THIS >
            static void call( T_THIS *_this, typename T_THIS :: t_Populator& _populator );
        };
        template<class T_THIS> void branch< Disc::unary_index > 
          :: call( T_THIS *_this, typename T_THIS :: t_Populator& _populator )
          {
            typedef typename Modifier::if_then_else
                    < 
                      boost::is_const<T_THIS> :: value, 
                      const typename T_THIS :: t_Derived,
                      typename T_THIS :: t_Derived
                    > :: type t_Derived;
            typedef typename boost::remove_const< t_Derived > :: type  t_NonConst;
            t_Derived *derived = static_cast< t_NonConst* >( _this  );
            t_Individual &indiv( *_populator );
            if( (*derived)( indiv ) ) indiv.invalidate();
          }
        template<> struct branch< Disc::const_unary_index > 
                      : public branch< Disc::unary_index > {};


        template<> struct branch< Disc::binary_index > 
        {
          template< class T_THIS >
            static void call( T_THIS *_this, typename T_THIS :: t_Populator& _populator );
        };
        template<class T_THIS> void branch< Disc::binary_index > 
          :: call( T_THIS *_this, typename T_THIS :: t_Populator& _populator )
          {
            typedef typename Modifier::if_then_else
                    < 
                      boost::is_const<T_THIS> :: value, 
                      const typename T_THIS :: t_Derived,
                      typename T_THIS :: t_Derived
                    > :: type t_Derived;
            typedef typename boost::remove_const< t_Derived > :: type  t_NonConst;
            t_Derived *derived = static_cast< t_NonConst* >( _this  );
            t_Individual &indiv1( *_populator );
            const t_Individual &indiv2( _populator.select() );
            if( (*derived)( indiv1, indiv2 ) ) indiv1.invalidate();
          }
        template<> struct branch< Disc::const_binary_index > 
                      : public branch< Disc::binary_index > {};

        template<> struct branch< Disc::other_index > 
        {
          template< class T_THIS >
            static void call( T_THIS *_this, typename T_THIS :: t_Populator& _populator );
        };
        template<class T_THIS> void branch< Disc::other_index > 
          :: call( T_THIS *_this, typename T_THIS :: t_Populator& _populator )
          {
            typedef typename Modifier::if_then_else
                    < 
                      boost::is_const<T_THIS> :: value, 
                      const typename T_THIS :: t_Derived,
                      typename T_THIS :: t_Derived
                    > :: type t_Derived;
            typedef typename boost::remove_const< t_Derived > :: type  t_NonConst;
            t_Derived *derived = static_cast< t_NonConst* >( _this  );
            (*derived)( _populator );
          }


      }
    }

  }
}



#endif 
