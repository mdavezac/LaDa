//
//  Version: $Id$
//
#ifndef _LADA_GA_PURELAYERS_MUTATION_H_
#define _LADA_GA_PURELAYERS_MUTATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/bind.hpp>

#include <opt/types.h>
#include <print/xmg.h>

#include "../bitstring/mutation.h"

namespace LaDa
{
  namespace GA
  {
    namespace PureLayers
    {
      //! Performs a bitstring crossover.
      template< class T_INDIVIDUAL, class T_CONCENTRATION >
        bool mutation( T_INDIVIDUAL& _o, types::t_real &_rate,
                       const T_CONCENTRATION &_concentration )
        {
          bool has_changed = GA::BitString::mutation< T_INDIVIDUAL >( _o, _rate );
          _concentration( _o.Object() );
          return has_changed;
        }

      namespace Factory
      {
        //! Creates a crossover operation.
        template< class T_FACTORY, class T_CONCENTRATION >
          void mutation( T_FACTORY &_factory,
                          boost::function<void( typename T_FACTORY::t_Populator& )>& _function,
                          const TiXmlElement &_node, const T_CONCENTRATION &_concentration )
          {
            typedef typename T_FACTORY :: t_Individual t_Individual;
            typedef typename T_FACTORY :: t_Populator t_Populator;
            typedef typename t_Individual :: t_IndivTraits :: t_Object t_Object;
            types::t_real rate(5e-1);
            if( _node.Attribute( "rate" ) )
              rate = boost::lexical_cast<types::t_real>( _node.Attribute("rate") );
            if( Fuzzy::leq( rate, 0e0 ) or Fuzzy::geq( rate, 1e0 ) ) rate = 5e-1;
            Operator::MakePopulator<t_Individual, t_Populator>::transform_unary
            ( 
              boost::bind( &PureLayers::mutation<t_Individual>, _1, rate,
                           boost::cref( _concentration ) ),
              _function 
            );
            Print::xmg << Print::Xmg::comment << "Mutation with rate=" << rate << Print::endl;
          }
      }

      namespace AddressOf
      {
        //! returns address of crossover operation.
        template< class T_FACTORY, class T_CONCENTRATION >
          void ( *mutation( const T_FACTORY&, const T_CONCENTRATION& ) )
                 ( T_FACTORY &, boost::function<void( typename T_FACTORY::t_Populator& )>&,
                   const TiXmlElement &, const T_CONCENTRATION& )
          { return &GA::PureLayers::Factory::mutation<T_FACTORY, T_CONCENTRATION>; }
      }
    } // namespace AlloyLayers
  } // namespace GA
} // namespace LaDa
#endif
