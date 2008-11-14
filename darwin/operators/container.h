//
//  Version: $Id$
//
#ifndef _LADA_DARWIN_OPERATORS_CONTAINER_H_
#define _LADA_DARWIN_OPERATORS_CONTAINER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <string>
#include <boost/function.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

#include <eo/utils/eoRNG.h>
#include <tinyxml/tinyxml.h>

#include <print/xmg.h>
#include <opt/types.h>

#include "make_populator.h"
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
      //! Creates a container which sequentially calls all of its functors.
      template< class T_INDIVIDUAL,
                class T_POPULATOR = eoPopulator< T_INDIVIDUAL > >
        class BaseContainer
        {
            //! Type of the functor container.
            typedef std::vector< boost::function<void(T_POPULATOR&) > > t_Functors;
            //! Small class to make Discriminate simpler.
            template< class T_FUNCTOR >
              class isother : public IsOther< T_INDIVIDUAL, T_FUNCTOR > {};
            //! Tag for discriminating between functors.
            template< bool D > class TypeTag {};
          public:
            //! Type of the indidividual
            typedef T_INDIVIDUAL t_Individual;
            //! Type of the populator
            typedef T_POPULATOR t_Populator;
     
            //! Constructor.
            BaseContainer() : functors_( new t_Functors ) {}
            //! Copy Constructor.
            BaseContainer( const BaseContainer& _c ) : functors_( _c.functors_ ) {}
            //! Virtual destrcutor, just to make sure.
            virtual ~BaseContainer() {}
     
            //! Functor over populator. Branches to correct format for t_Derived.
            void operator()( t_Populator& _populator )
            {
              typename t_Functors :: iterator i_functor = functors_->begin();
              typename t_Functors :: iterator i_functor_end = functors_->end();
              for(; i_functor != i_functor_end; ++i_functor )
                (*i_functor)( _populator );
            }
            //! Connects the callback.
            template< class T_FUNCTOR >
              void connect( T_FUNCTOR& _functor )
              { 
                typename t_Functors :: value_type popfunc;
                MakePopulator< t_Individual, t_Populator > :: transform( _functor, popfunc );
                functors_->push_back( typename t_Functors :: value_type( popfunc ) );
              }

            //! Connects the callback.
            template< class T_FUNCTOR >
              void connect( boost::shared_ptr<T_FUNCTOR>& _functor )
              {
                typename t_Functors :: value_type popfunc;
                MakePopulator< t_Individual, t_Populator > :: transform( *_functor, popfunc );
                functors_->push_back( popfunc );
              }

          protected:
            //! The functors.
            boost::shared_ptr<t_Functors> functors_;
            
        };

      //! Creates a container class with some selection ability over which functors to call.
      template< class T_INDIVIDUAL,
                class T_OPSELECTIONPOLICY,
                class T_POPULATOR = eoPopulator< T_INDIVIDUAL > >
        class Container : protected BaseContainer< T_INDIVIDUAL, T_POPULATOR >
        {
            //! Type of the functor container.
            typedef std::vector< boost::function<void(T_POPULATOR&) > > t_Functors;
            //! Type of the rates container.
            typedef std::vector< types::t_real > t_Rates;
            //! type of the base class;
            typedef BaseContainer<T_INDIVIDUAL, T_POPULATOR>  t_Base;
          public:
            //! Type of the indidividual
            typedef T_INDIVIDUAL t_Individual;
            //! Type of the operator selection policy.
            typedef T_OPSELECTIONPOLICY t_OpSelectionPolicy;
            //! Type of the populator
            typedef T_POPULATOR t_Populator;
     
            //! Constructor.
            Container() : t_Base(), rates_(new t_Rates) {}
            //! Copy Constructor.
            Container( const Container& _c ) : t_Base(_c), rates_(_c.rates_) {}
            //! Virtual destrcutor, just to make sure.
            virtual ~Container() {}
     
            //! Functor over populator. Branches to correct format for t_Derived.
            void operator()( t_Populator& _populator )
            {
              __ASSERT( not rates_.get(), "rates_ container does not exist.\n" )
              __ASSERT( not functors_.get(), "functors_ container does not exist.\n" )
              __ASSERT( rates_->size() != functors_->size(),
                        "functors_ and rates_ containers have different sizes.\n" )
              t_OpSelectionPolicy selection( *rates_ );
              typename t_Functors :: iterator i_functor = functors_->begin();
              typename t_Functors :: iterator i_functor_end = functors_->end();
              t_Rates :: const_iterator i_rate = rates_->begin();
              for(; i_functor != i_functor_end; ++i_functor, ++i_rate )
                if( selection( *i_rate ) ) (*i_functor)( _populator );
            }
            //! Connects the callback.
            template< class T_FUNCTOR >
              void connect( T_FUNCTOR& _functor, types::t_real _rate )
              {
                connect( _functor ); 
                rates_->push_back( _rate );
              }
            //! Connects the callback.
            template< class T_FUNCTOR >
              void connect( boost::shared_ptr<T_FUNCTOR>& _functor,
                            types::t_real _rate )
              {
                connect( *_functor ); 
                rates_->push_back( _rate );
              }
            
          private:
            using t_Base :: functors_;
            using t_Base :: connect;
            //! The rates.
            boost::shared_ptr< t_Rates > rates_;
        };

      //! A policy for a sequential container of operators.
      struct Sequential
      {
        // Constructor.
        Sequential( const std::vector< types::t_real >& ) {}
        //! Functor.
        bool operator()( types::t_real _rate ) { return eo::rng.flip( _rate ); }
        //! Name of the container/policy.
        static const std::string name;
      };

      //! A policy for a proportional container of operators.
      struct Proportional
      {
        // Constructor.
        Proportional   ( const std::vector< types::t_real >&  _rates )
                     : i( eo::rng.roulette_wheel(_rates) ), index(0) {}
        //! Functor.
        bool operator()( types::t_real _rate ) { return index++ == i; }

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
                        const TiXmlElement &_node )
        {
          typedef typename T_FACTORY :: t_Individual t_Individual;
          typedef typename T_FACTORY :: t_Populator t_Populator;
          typedef Operator::Container<t_Individual, T_POLICY, t_Populator> t_Container;
          t_Container container;
          const TiXmlElement *child = _node.FirstChildElement();
          Print::xmg << Print::Xmg::comment << T_POLICY::name << Print::endl
                     << Print::Xmg::indent;
          size_t nbcreated(0);
          for(; child; child = child->NextSiblingElement() )
          {
            const std::string key( child->Value() );
            if( not _factory.exists( key ) ) continue;
            types::t_real prob(1e0);
            if( child->Attribute( "prob" ) )
              prob = boost::lexical_cast<types::t_real>( child->Attribute( "prob" ) );
            boost::function<void(t_Populator&)> function;
            _factory( key, function, *child );
            Print::xmg << Print::Xmg::addtolast << " prob = " << prob << Print::endl;
            container.connect( function, prob );
            ++nbcreated;
          }
          if( nbcreated == 0 )
          {
            Print::xmg << Print::Xmg::removelast;
            return;
          }
          Print::xmg << Print::Xmg :: unindent;
          _function = container;
        }

      //! Specializes void container() to a sequential container.
      template< class T_FACTORY >
        void sequential( T_FACTORY &_factory,
                         boost::function<void( typename T_FACTORY::t_Populator& )>&
                           _function,
                         const TiXmlElement &_node )
        {
          container< T_FACTORY, Operator::Sequential>( _factory, _function, _node );
        }
      //! Specializes void container() to a proportional container.
      template< class T_FACTORY >
        void proportional( T_FACTORY &_factory,
                            boost::function<void( typename T_FACTORY::t_Populator& )>&
                              _function,
                            const TiXmlElement &_node )
        {
          container<T_FACTORY, Operator::Proportional>( _factory, _function, _node );
        }

      template< class T_FACTORY >
        void containers( T_FACTORY &_factory,
                         boost::function<void( typename T_FACTORY::t_Populator& )>&
                           _function,
                         const TiXmlElement &_node )
        {
          std::string type("and");
          if( _node.Attribute("type") ) type = _node.Attribute("type");
          __DOASSERT( type != "and" and type != "or",
                      "Tag Operators expects an attribute \"type\" "
                      "with values \"and\" or \"or\".\n" )
          type == "and" ?
            container<T_FACTORY, Operator::Sequential>( _factory, _function, _node ):
            container<T_FACTORY, Operator::Proportional>( _factory, _function, _node );
        }

    } // Factory namespace 
  } // GA namespace 
} // LaDa namespace 



#endif 
