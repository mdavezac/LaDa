//
//  Version: $Id$
//
#ifndef _LADA_GA_ALLOYLAYERS_CROSSOVER_H_
#define _LADA_GA_ALLOYLAYERS_CROSSOVER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../bitstring/crossover.h"

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
          _o.random_symmetry();
          return BitString::crossover( _o, _p, _rate );
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
          if( Fuzzy::leq( rate, 0e0 ) or Fuzzy::geq( rate, 1e0 ) ) rate = 5e-1;
          Operator::MakePopulator<t_Individual, t_Populator>::transform_binary
          ( 
            boost::bind( &crossover<t_Individual>, _1, _2, rate ),
            _function 
          );
          Print::xmg << Print::Xmg::comment << "Crossover with rate=" << rate << Print::endl;
        }
    } // namespace AlloyLayers
  } // namespace GA
} // namespace LaDa
#endif
