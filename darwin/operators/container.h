//
//  Version: $Id$
//
#ifndef _LADA_DARWIN_OPERATORS_CONTAINER_H_
#define _LADA_DARWIN_OPERATORS_CONTAINER_H_

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
      //! Creates a unary callback operator.
      template< class T_INDIVIDUAL,
                class T_OPSELECTIONPOLICY,
                class T_POPULATOR = eoPopulator< T_INDIVIDUAL > >
        class Container
        {
            //! Type of the functor container.
            typedef std::vector< boost::function<void(t_Populator&) > > t_Functors;
            //! Type of the rates container.
            typedef std::vector< types::t_real > t_Rates;
          public:
            //! Type of the indidividual
            typedef T_INDIVIDUAL t_Individual;
            //! Type of the operator selection policy.
            typedef T_OPSELECTIONPOLICY t_OpSelectionPolicy;
            //! Type of the populator
            typedef T_POPULATOR t_Populator;
     
            //! Virtual destrcutor, just to make sure.
            virtual ~Container() {}
     
            //! Functor over populator. Branches to correct format for t_Derived.
            void operator( t_Populator& _populator )
            {
              t_OpSelectionPolicy selection( rates_ );
              t_Functors :: iterator i_functor = functors_.begin();
              t_Functors :: iterator i_functor_end = functors_.end();
              for(; i_functor != i_functor_end; ++i_functor, ++i_rate )
                if( selection( *i_rate ) ) (*i_functor)( _populator );
            }
            //! Connects the callback.
            template< class T_FUNCTOR > connect( T_FUNCTOR& _functor, types::t_real _rate )
              { connect( _functor, _rate, TypeTag< Discriminate<t_Individual>::value >() ); }
            //! Connects the callback.
            template< class T_FUNCTOR > connect( boost::shared_ptr<T_FUNCTOR>& _functor,
                                                 types::t_real _rate )
              { connect( _functor, _rate ); }
            
          private:
            //! Branch for unary operators.
            template< class T_FUNCTOR >
              connect( T_FUNCTOR& _functor, types::t_real _rate, TypeTag<2>& )
              {
                PopulatorOperator< T_FUNCTOR, t_Individual, t_Populator > popfunc( _functor );
                connect( popfunc, types::t_real, TypeTag<1> );
              }
            //! Branch for binary operators.
            template< class T_FUNCTOR >
              connect( T_FUNCTOR& _functor, types::t_real _rate, TypeTag<3>& );
              {
                PopulatorOperator< T_FUNCTOR, t_Individual, t_Populator > popfunc( _functor );
                connect( popfunc, types::t_real, TypeTag<1> );
              }
            //! Branch for populator operator (asssumed, not detected).
            template< class T_FUNCTOR >
              connect( T_FUNCTOR& _functor, types::t_real _rate, TypeTag<1>& )
              {
                functors_.push_back( t_Functors :: value_type( _functor ) );
                rates_.push_back( _rate );
              }

            //! The functors.
            t_Functors functors_;
            //! The rates.
            t_Rates rates_;
        };

      //! A policy for a sequential container of operators.
      struct Sequential
      {
        // Constructor.
        Sequential( const std::vector< types::t_real >& ) {}
        //! Functor.
        bool operator()( types::t_real& _rate ) { eo::rng.flip( _rate ); }
      };

      //! A policy for a proportional container of operators.
      struct Proportional
      {
        // Constructor.
        Sequential   ( const std::vector< types::t_real >&  _rates )
                   : i( _rates ), index(0) {}
        //! Functor.
        bool operator()( types::t_real& _rate ) { return index++ == i; }

        //! The value for which to accept the operation.
        types::t_unsigned i;
        //! The current index (or number of calls).
        types::t_unsigned index;
      };

    } // operator namespace.

    namespace Factory
    {
      template< class T_FACTORY, class T_INDIVIDUAL, class T_POLICY class T_POPULATOR >
        boost::shared_ptr< Operator::Container< T_INDIVIDUAL, T_POLICY, T_POPULATOR > >
          container( T_FACTORY &_factory, TiXmlElement &_node )
          {
            typedef Operator::Container< T_INDIVIDUAL, T_POLICY, T_POPULATOR > t_Container;
            typedef boost::shared_ptr< t_Container > t_Result;
            t_Result result( new t_Container );
            const TiXmlElement *child = _node.FirstChildElement();
            for(; child; child = child->NextSiblingElement() )
            {
              if( not _factory.exists( child ) ) continue;
              types::t_real rate(1e0);
              if( child->Attribute( "rate" ) )
                rate = boost::lexical_cast<types::t_real>( child->Attribute( "rate" ) );
              result->connect( _factory( _child ), rate );
            }
            return result;
          }

      template< class T_FACTORY, class T_INDIVIDUAL, class T_POPULATOR >
        boost::shared_ptr< Operator::Container< T_INDIVIDUAL, Operator::Sequential, T_POPULATOR > >
          sequential( T_FACTORY &_factory, TiXmlElement &_node )
           { return container< T_FACTORY, T_INDIVIDUAL, Operator::Sequential, T_POPULATOR >
                             ( _factory, _node ); }
      template< class T_FACTORY, class T_INDIVIDUAL, class T_POPULATOR >
        boost::shared_ptr< Operator::Container< T_INDIVIDUAL, Operator::Proportional, T_POPULATOR > >
          sequential( T_FACTORY &_factory, TiXmlElement &_node )
           { return container< T_FACTORY, T_INDIVIDUAL, Operator::Proportional, T_POPULATOR >
                             ( _factory, _node ); }

    } // Factory namespace 
  } // GA namespace 
} // LaDa namespace 



#endif 
