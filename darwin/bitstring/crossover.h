//
//  Version: $Id$
//
#ifndef _LADA_GA_BITSTRING_CROSSOVER_H_
#define _LADA_GA_BITSTRING_CROSSOVER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/bind.hpp>

#include <opt/types.h>
#include <print/xmg.h>

#include "../operators/make_populator.h"

namespace LaDa
{
  namespace GA
  {
    namespace BitString
    {
      //! Performs a bitstring crossover.
      template< class T_INDIVIDUAL >
        bool crossover( T_INDIVIDUAL& _o, const T_INDIVIDUAL &_p, types::t_real &_rate )
        {
          typedef typename T_INDIVIDUAL :: t_IndivTraits :: t_Object t_Object;
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

      template< class T_FACTORY >
        void crossover_factory( T_FACTORY &_factory,
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
          if( math::leq( rate, 0e0 ) or math::geq( rate, 1e0 ) ) rate = 5e-1;
          Operator::MakePopulator<t_Individual, t_Populator>::transform_binary
          ( 
            boost::bind( &crossover<t_Individual>, _1, _2, rate ),
            _function 
          );
          Print::xmg << Print::Xmg::comment << "Crossover with rate=" << rate << Print::endl;
        }
    } // namespace BitString
  } // namespace GA
} // namespace LaDa
#endif
