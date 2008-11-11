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
            template< class T_FUNCTOR >
              void connect( T_FUNCTOR& _functor, types::t_real _rate )
              { connect( _functor, _rate,
                         TypeTag< Discriminate<t_Individual>::value >() ); }
            //! Connects the callback.
            template< class T_FUNCTOR >
              void connect( boost::shared_ptr<T_FUNCTOR>& _functor,
                            types::t_real _rate )
              { connect( *_functor, _rate ); }
            
          private:
            //! Branch for unary operators.
            template< class T_FUNCTOR >
              connect( T_FUNCTOR& _functor, types::t_real _rate, TypeTag<2>& )
              {
                PopulatorOperator< T_FUNCTOR, t_Individual,
                                   t_Populator > popfunc( _functor );
                connect( popfunc, types::t_real, TypeTag<1> );
              }
            //! Branch for binary operators.
            template< class T_FUNCTOR >
              connect( T_FUNCTOR& _functor, types::t_real _rate, TypeTag<3>& );
              {
                PopulatorOperator< T_FUNCTOR, t_Individual,
                                   t_Populator > popfunc( _functor );
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
        //! Name of the container/policy.
        static const std::string name;
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
        //! Name of the container/policy.
        static const std::string name;
      };

    } // operator namespace.

    namespace Factory
    {
      template< class T_FACTORY,class T_POLICY >
        void container( T_FACTORY &_factory,
                        boost::function<void(typename T_FACTORY::t_Populator&)>&
                          _function,
                        TiXmlElement &_node )
        {
          typedef Operator::Container
                  < 
                    typename T_FACTORY :: t_Individual, 
                    T_POLICY, 
                    typename T_FACTORY :: t_Populator 
                  > t_Container;
          typedef boost::shared_ptr< t_Container > t_Result;
          t_Result result( new t_Container );
          const TiXmlElement *child = _node.FirstChildElement();
          Print::xmg << Print::Xmg::comment << T_POLICY::name << Print::endl
                     << Print::Xmg::indent;
          for(; child; child = child->NextSiblingElement() )
          {
            if( not _factory.exists( child->Value() ) ) continue;
            types::t_real rate(1e0);
            if( child->Attribute( "prob" ) )
              rate = boost::lexical_cast<types::t_real>( child->Attribute( "prob" ) );
            boost::function<void(T_POPULATOR&)> function;
            _factory( function, *child );
            Print::xmg << Print::Xmg::addtolast << " prob = " << prob << Print::endl;
            result->connect( function, rate );
          }
          Print::xmg << Print::Xmg :: unindent;
          return result;
        }

      //! Specializes void container() to a sequential container.
      template< class T_FACTORY >
        void sequential( T_FACTORY &_factory,
                         boost::function<void( typename T_FACTORY::t_Populator& )>&
                           _function,
                         const TiXmlElement &_node )
        {
          container
          < 
            T_FACTORY, 
            Operator::Sequential,
            typename T_FACTORY :: t_Populator
          > ( _factory, _function, _node );
        }
      //! Specializes void container() to a proportional container.
      template< class T_FACTORY >
        void proportional( T_FACTORY &_factory,
                            boost::function<void( typename T_FACTORY::t_Populator& )>&
                              _function,
                            const TiXmlElement &_node )
        {
          container
          < 
            T_FACTORY, 
            Operator::Sequential
          > ( _factory, _function, _node );
        }
    } // Factory namespace 
  } // GA namespace 
} // LaDa namespace 



#endif 
