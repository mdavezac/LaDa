//
//  Version: $Id$
//
#ifndef _LADA_GA_PURELAYERS_CROSSOVER_H_
#define _LADA_GA_PURELAYERS_CROSSOVER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/bind.hpp>

#include <opt/types.h>
#include <print/xmg.h>

#include "../bitstring/crossover.h"

namespace LaDa
{
  namespace GA
  {
    namespace PureLayers
    {
      //! Performs a bitstring crossover.
      template< class T_INDIVIDUAL, class T_CONCENTRATION >
        bool crossover( T_INDIVIDUAL& _o, const T_INDIVIDUAL &_p, types::t_real &_rate,
                        const T_CONCENTRATION &_concentration )
        {
          bool has_changed = GA::BitString::crossover< T_INDIVIDUAL >( _o, _p, _rate );
          _concentration( _o.Object() );
          return has_changed;
        }

      namespace Factory
      {
        //! Creates a crossover operation.
        template< class T_FACTORY, class T_CONCENTRATION >
          void crossover( T_FACTORY &_factory,
                          boost::function<void( typename T_FACTORY::t_Populator& )>& _function,
                          const TiXmlElement &_node, const T_CONCENTRATION &_concentration )
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
              boost::bind( &PureLayers::crossover<t_Individual>, _1, _2,
                           rate, boost::cref( _concentration ) ),
              _function 
            );
            Print::xmg << Print::Xmg::comment << "Crossover with rate=" << rate << Print::endl;
          }
      }

      namespace AddressOf
      {
        //! returns address of crossover operation.
        template< class T_FACTORY, class T_CONCENTRATION >
          void ( *crossover( const T_FACTORY&, const T_CONCENTRATION& ) )
                 ( T_FACTORY &, boost::function<void( typename T_FACTORY::t_Populator& )>&,
                   const TiXmlElement &, const T_CONCENTRATION& )
          { return &GA::PureLayers::Factory::crossover<T_FACTORY, T_CONCENTRATION>; }
      }
    } // namespace AlloyLayers
  } // namespace GA
} // namespace LaDa
#endif
