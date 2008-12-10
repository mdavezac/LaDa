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

#include "../operators/make_populator.h"

namespace LaDa
{
  namespace GA
  {
    namespace AlloyLayers
    {
      namespace Factory
      {
        template< class T_FACTORY, class T_EVALUATOR >
          void random( T_FACTORY &_factory,
                       boost::function<void( typename T_FACTORY::t_Populator& )>&
                         _function,
                       const TiXmlElement &_node, 
                       T_EVALUATOR &_evaluator )
          {
            typedef typename T_FACTORY :: t_Individual t_Individual;
            typedef typename T_FACTORY :: t_Populator t_Populator;
            Operator::MakePopulator<t_Individual, t_Populator>::transform_unary
            ( 
              boost::bind( &T_EVALUATOR::initialize, boost::ref(_evaluator), _1 ),
              _function 
            );
            Print::xmg << Print::Xmg::comment << "Random" << Print::endl;
          }
      } // namespace Factory

      namespace AddressOf
      {
        //! Helper function to return address of GA::AlloyLayers::Factory::random().
        template< class T_FACTORY, class T_EVALUATOR >
          void(* random( const T_FACTORY&, const T_EVALUATOR& ) )
                  ( T_FACTORY &, boost::function<void( typename T_FACTORY::t_Populator& )>&,
                    const TiXmlElement &, T_EVALUATOR & )
          { return &GA::AlloyLayers::Factory::random< T_FACTORY, T_EVALUATOR >; }
      }
     
    } // namespace AlloyLayers
  } // namespace GA
} // namespace LaDa
#endif
