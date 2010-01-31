//
//  Version: $Id$
//
#ifndef _LADA_GA_ALLOYLAYERS_MUTATION_H_
#define _LADA_GA_ALLOYLAYERS_MUTATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/bind.hpp>
#include<boost/ref.hpp>

#include <opt/types.h>
#include <print/xmg.h>

#include "../operators/make_populator.h"

namespace LaDa
{
  namespace GA
  {
    namespace AlloyLayers
    {
      //! Performs a bitstring mutation.
      template< class T_INDIVIDUAL, class T_EVALUATOR >
        bool mutation( T_INDIVIDUAL& _o, types::t_real &_rate, T_EVALUATOR const &_evaluator )
        {
          typedef typename T_INDIVIDUAL :: t_IndivTraits :: t_Object t_Object;
          typedef typename t_Object :: t_Container :: value_type t_Type;
          bool has_changed = false;
          Crystal::Structure structure( _evaluator.get_structure() );
          _evaluator.translate( _o, structure );
          foreach( Crystal::Structure::t_Atom atom, structure.atoms )
          {
            if( not rng.flip( _rate ) ) continue;
            atom.type -= atom.type;
            has_changed = true;
          }
          _evaluator.translate(structure, _o);
          return has_changed;
        }

      //! Factory functor for mutations.
      template< class T_EVALUATOR >
        struct MutationFactory
        {
           MutationFactory( T_EVALUATOR const &_eval ) : evaluator(_eval) {}
           MutationFactory( MutationFactory const &_c ) : evaluator(_c.evaluator) {}
           T_EVALUATOR const &evaluator;
           typedef void result_type;
           template< class T_FACTORY>
             void operator()( T_FACTORY &_factory,
                              boost::function<void( typename T_FACTORY::t_Populator& )>&
                                _function,
                              const TiXmlElement &_node ) const
             {
               typedef typename T_FACTORY :: t_Individual t_Individual;
               typedef typename T_FACTORY :: t_Populator t_Populator;
               typedef typename t_Individual :: t_IndivTraits :: t_Object t_Object;
               types::t_real rate(5e-1);
               if( _node.Attribute( "rate" ) )
                 rate = boost::lexical_cast<types::t_real>( _node.Attribute("rate") );
               if( math::leq( rate, 0e0 ) or math::geq( rate, 1e0 ) ) rate = 5e-1;
               Operator::MakePopulator<t_Individual, t_Populator>::transform_unary
               ( 
                 boost::bind(  &mutation<t_Individual, T_EVALUATOR>, _1, rate,
                               boost::cref(evaluator) ),
                 _function 
               );
               Print::xmg << Print::Xmg::comment << "Mutation with rate=" << rate << Print::endl;
             }
        };
      //! Creates mutation factory object.
      template< class T_EVALUATOR >
        MutationFactory<T_EVALUATOR> mutation_factory( T_EVALUATOR const &_evaluator )
          { return MutationFactory<T_EVALUATOR>( _evaluator ); }

    } // namespace AlloyLayers
  } // namespace GA
} // namespace LaDa
#endif
