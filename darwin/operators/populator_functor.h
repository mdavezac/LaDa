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
      template< class T_DERIVED, class T_INDIVIDUAL, class T_POPULATOR > class PopulatorFunctor;
      //! \endcond

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
            //! A tag structure.
            template< size_t  D > struct TypeTag {};
            //! Type of the discrimination metafunction.
            typedef Discriminate<T_INDIVIDUAL, T_DERIVED> t_Discriminate;
            //! "Discrimination" index of the functor.
            const static size_t discrimination = t_Discriminate :: value;
     
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
              { branch( _populator, TypeTag< discrimination >() ); }
            //! Functor over populator. Branches to correct format for t_Derived.
            void operator()( t_Populator& _populator ) const
              { branch( _populator, TypeTag< discrimination >() ); }
            
          private:
            //! Type of the tag for "other" functors.
            typedef TypeTag< t_Discriminate :: other_index > other_tag;
            //! Type of the tag for unary functors.
            typedef TypeTag< t_Discriminate :: unary_index > unary_tag;
            //! Type of the tag for binary functors.
            typedef TypeTag< t_Discriminate :: binary_index > binary_tag;
            //! Type of the tag for const unary functors.
            typedef TypeTag< t_Discriminate :: const_unary_index > const_unary_tag;
            //! Type of the tag for const binary functors.
            typedef TypeTag< t_Discriminate :: const_binary_index > const_binary_tag;
            //! Unary overload.
            void branch( t_Populator& _populator, const unary_tag& )
            {
              t_Derived *_this = static_cast< t_Derived* >( this  );
              t_Individual &indiv( *_populator );
              if( (*_this)( indiv ) ) indiv.invalidate();
            }
            //! Binary overload.
            void branch( t_Populator& _populator, const binary_tag& )
            {
              t_Derived *_this = static_cast< t_Derived* >( this  );
              t_Individual &indiv1( *_populator );
              const t_Individual &indiv2( _populator.select() );
              if( (*_this)( indiv1, indiv2 ) ) indiv1.invalidate();
            }
            //! "Other" overload.
            void branch( t_Populator& _populator, const other_tag& )
              { t_Derived :: operator()( _populator ); }
            //! const Unary overload.
            void branch( t_Populator& _populator, const const_unary_tag& ) const
            {
              const t_Derived *_this = static_cast< const t_Derived* >( this  );
              t_Individual &indiv( *_populator );
              if( (*_this)( indiv ) ) indiv.invalidate();
            }
            //! const Binary overload.
            void branch( t_Populator& _populator, const const_binary_tag& ) const
            {
              const t_Derived *_this = static_cast< const t_Derived* >( this  );
              t_Individual &indiv1( *_populator );
              const t_Individual &indiv2( _populator.select() );
              if( (*_this)( indiv1, indiv2 ) ) indiv1.invalidate();
            }
            //! const "Other" overload.
            void branch( t_Populator& _populator, const other_tag& ) const
              { t_Derived :: operator()( _populator ); }
            //! Assigns this object to \a _function.
            void assign( boost::function< void( t_Populator& ) >& );
        };

      template< class T_DERIVED, class T_INDIVIDUAL, class T_POPULATOR >
        void assign( PopulatorFunctor<T_DERIVED, T_INDIVIDUAL, T_POPULATOR>& _b, 
                     boost::function< void( T_POPULATOR& ) > &_a )
        {
          typedef PopulatorFunctor< T_DERIVED, T_INDIVIDUAL, T_POPULATOR > t_PopFunc;
          const T_DERIVED * const derived = static_cast< const T_DERIVED* const >( &_b );
          void (t_PopFunc::*ptr_popop)( T_POPULATOR& ) = &t_PopFunc :: operator();
          _a = boost::bind( ptr_popop, *derived, _1 );
        }

      template< class T_DERIVED, class T_INDIVIDUAL, class T_POPULATOR >
        void assign( const PopulatorFunctor<T_DERIVED, T_INDIVIDUAL, T_POPULATOR>& _b, 
                     boost::function< void( T_POPULATOR& ) > &_a )
        {
          typedef PopulatorFunctor< T_DERIVED, T_INDIVIDUAL, T_POPULATOR > t_PopFunc;
          const T_DERIVED * const derived = static_cast< const T_DERIVED* const >( &_b );
          void (t_PopFunc::*ptr_popop)( T_POPULATOR& ) const = &t_PopFunc :: operator();
          _a = boost::bind( ptr_popop, *derived, _1 );
        }

    }

  }
}



#endif 
