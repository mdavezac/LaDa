//
//  Version: $Id$
//
#ifndef _LADA_GA_ALLOYLAYERS_OPERATORS_H_
#define _LADA_GA_ALLOYLAYERS_OPERATORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/bind.hpp>

#include <opt/types.h>
#include <print/xmg.h>

#include "../operators/callback.h"

namespace LaDa
{
  namespace GA
  {
    namespace AlloyLayers
    {
      //! Performs a bitstring crossover.
      template< class T_INDIVIDUAL >
        bool crossover( T_INDIVIDUAL& _o, const T_INDIVIDUAL &_p, types::t_real &_rate )
        {
          typedef typename T_INDIVIDUAL :: t_Object t_Object;
          typename t_Object :: iterator  i_var_off = _o.Object().Container().begin();
          typename t_Object :: iterator  i_var_off_end = _o.Object().Container().end();
          typename t_Object :: const_iterator  i_var_par = _p.Object().Container().begin();
          bool has_changed = false;
          for(; i_var_off != i_var_off_end; ++i_var_off, ++i_var_par )
          {
            if( not rng.flip( _rate ) ) continue;
            *i_var_off = *i_var_par;
            has_changed = true;
          }
          return has_changed;
        }
      //! Performs a bitstring mutation.
      template< class T_INDIVIDUAL >
        bool mutation( T_INDIVIDUAL& _o, types::t_real &_rate )
        {
          typedef typename T_INDIVIDUAL :: t_Object t_Object;
          typedef typename t_Object :: t_Container :: value_type t_Type;
          typename t_Object :: iterator  i_var_off = _o.Object().Container().begin();
          typename t_Object :: iterator  i_var_off_end = _o.Object().Container().end();
          bool has_changed = false;
          for(; i_var_off != i_var_off_end; ++i_var_off )
          {
            if( not rng.flip( _rate ) ) continue;
            *i_var_off = t_Type(  rng.uniform( 2e0 ) - 1e0 ); 
            has_changed = true;
          }
          return has_changed;
        }

      namespace Factory
      {
        template< class T_FACTORY >
          void crossover( T_FACTORY &_factory,
                          boost::function<void( typename T_FACTORY::t_Populator& )>&
                            _function,
                          const TiXmlElement &_node )
          {
            typedef typename T_FACTORY :: t_Individual t_Individual;
            typedef typename T_FACTORY :: t_Populator t_Populator;
            typedef typename t_Individual :: t_IndivTraits :: t_Object t_Object;
            types::t_real rate(5e-1);
            if( _node.Attribute( "rate" ) )
              rate = boost::lexical_cast<types::t_real>( _node.Attribute("rate") );
            if( Fuzzy::leq( rate, 0e0 ) or Fuzzy::geq( rate, 1e0 ) ) rate = 5e-1;
            Operator::BinaryCallBack<t_Individual, t_Populator> callback;
            callback.connect
            (
              boost::bind( &::LaDa::GA::AlloyLayers::crossover<t_Object>,
                           _1, _2, rate )
            );
            Operator::assign( callback, _function );
            Print::xmg << Print::Xmg::comment << "Crossover with rate=" << rate << Print::endl;
          }
        template< class T_FACTORY >
          void mutation( T_FACTORY &_factory,
                         boost::function<void( typename T_FACTORY::t_Populator& )>&
                           _function,
                         const TiXmlElement &_node )
          {
            typedef typename T_FACTORY :: t_Individual t_Individual;
            typedef typename T_FACTORY :: t_Populator t_Populator;
            typedef typename t_Individual :: t_IndivTraits :: t_Object t_Object;
            types::t_real rate(5e-1);
            if( _node.Attribute( "rate" ) )
              rate = boost::lexical_cast<types::t_real>( _node.Attribute("rate") );
            if( Fuzzy::leq( rate, 0e0 ) or Fuzzy::geq( rate, 1e0 ) ) rate = 5e-1;
            Operator::BinaryCallBack<t_Individual, t_Populator> callback;
            callback.connect
            (
              boost::bind(  &::LaDa::GA::AlloyLayers::mutation<t_Object>, _1, rate )
            );
            Operator::assign( callback, _function );
            Print::xmg << Print::Xmg::comment << "Mutation with rate=" << rate << Print::endl;
          }

        template< class T_FACTORY, class T_EVALUATOR >
          void random( T_FACTORY &_factory,
                       boost::function<void( typename T_FACTORY::t_Populator& )>&
                         _function,
                       const TiXmlElement &_node, 
                       T_EVALUATOR &_evaluator )
          {
            typedef typename T_FACTORY :: t_Individual t_Individual;
            typedef typename T_FACTORY :: t_Populator t_Populator;
            Operator::UnaryCallBack<t_Individual, t_Populator> callback;
            callback.connect
            (
              boost::bind( &T_EVALUATOR::initialize, boost::ref(_evaluator), _1 )
            );
            Operator::assign( callback, _function );
            Print::xmg << Print::Xmg::comment << "Random" << Print::endl;
          }
     
      } // namespace Factory
    } // namespace AlloyLayers
  } // namespace GA
} // namespace LaDa
#endif
