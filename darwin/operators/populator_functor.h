//
//  Version: $Id$
//
#ifndef _LADA_DARWIN_OPERATORS_POPULATOR_FUNCTOR_H_
#define _LADA_DARWIN_OPERATORS_POPULATOR_FUNCTOR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "discriminate.h"

namespace LaDa
{
  namespace GA
  {
    //! Holds general stuff for GA operators.
    namespace Operator
    {
      //! Creates a genetic operator, eg, a functor of a T_POPULATOR
      template< class T_DERIVED,
                class T_INDIVIDUAL = typename T_DERIVED :: t_Individual,
                class T_POPULATOR = eoPopulator< T_INDIVIDUAL > >
        class PopulatorFunctor
        {
          private:
            //! A tag structure.
            template< int D > struct TypeTag {};
     
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
            void operator( t_Populator& _populator )
              { branch( _populator, TypeTag< Discriminate<t_Derived> :: value >() ); }
            
          private:
            //! Unary overload.
            void branch( t_Populator& _populator, TypeTag< 2 >& )
            {
              t_Derived *_this = static_cast< t_Derived* >( *this  );
              t_Individual &indiv( *_populator );
              if( (*_this)( indiv ) ) invid.invalidate();
            }
            //! Binary overload.
            void branch( t_Populator& _populator, TypeTag< 3 >& )
            {
              t_Derived *_this = static_cast< t_Derived* >( *this  );
              t_Individual &indiv1( *_populator );
              const t_Individual &indiv2( _populator.select() );
              if( (*_this)( indiv1, indiv2 ) ) invid.invalidate();
            }
        };

    }

  }
}



#endif 
