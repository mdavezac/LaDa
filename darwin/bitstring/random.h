//
//  Version: $Id$
//
#ifndef _LADA_GA_BISTRING_CROSSOVER_H_
#define _LADA_GA_BISTRING_CROSSOVER_H_

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
      //! Creates a random individual
      template< class T_OBJECT >
        bool random( T_OBJECT& _o,
                     typename T_OBJECT::t_Container::value_type min,
                     typename T_OBJECT::t_Container::value_type max )
        {
          typename t_Object :: iterator  i_var_off = _o.Container().begin();
          typename t_Object :: iterator  i_var_off_end = _o.Container().end();
          bool has_changed = false;
          for(; i_var_off != i_var_off_end; ++i_var_off )
            *i_var_off = rng.uniform( max - min ) - min;
          return true;
        }

      //! Creates a random two valued individual.
      template< class T_OBJECT >
        bool random_two_valued( T_OBJECT& _o,
                                typename T_OBJECT::t_Container::value_type _one,
                                typename T_OBJECT::t_Container::value_type _two )
        {
          typename t_Object :: iterator  i_var_off = _o.Container().begin();
          typename t_Object :: iterator  i_var_off_end = _o.Container().end();
          bool has_changed = false;
          for(; i_var_off != i_var_off_end; ++i_var_off )
            *i_var_off = rng.flip() ? _one: _two;
          return true;
        }

      template< class T_FACTORY >
        void random_factory( T_FACTORY &_factory,
                             boost::function<void( typename T_FACTORY::t_Populator& )>&
                               _function,
                             const TiXmlElement &_node,
                             typename T_FACTORY::t_Individual::t_IndivTraits
                                               ::t_Object::t_Type _max, 
                             typename T_FACTORY::t_Individual::t_IndivTraits
                                               ::t_Object::t_Type _min )
        {
          typedef typename T_FACTORY :: t_Individual t_Individual;
          typedef typename T_FACTORY :: t_IndivTraits :: t_Object t_Object;
          typedef typename T_FACTORY :: t_Populator t_Populator;
          Operator::MakePopulator<t_Individual, t_Populator>::transform_unary
          ( 
            boost::bind( &random<t_Object>, _1, _max, _min ),
            _function 
          );
          Print::xmg << Print::Xmg::comment << "Random" << Print::endl;
        }

      template< class T_FACTORY >
        void random_two_valued_factory( T_FACTORY &_factory,
                                        boost::function<void( typename T_FACTORY::t_Populator& )>&
                                          _function,
                                        const TiXmlElement &_node,
                                        typename T_FACTORY::t_Individual::t_IndivTraits
                                                          ::t_Object::t_Type _one, 
                                        typename T_FACTORY::t_Individual::t_IndivTraits
                                                          ::t_Object::t_Type _two )
        {
          typedef typename T_FACTORY :: t_Individual t_Individual;
          typedef typename T_FACTORY :: t_IndivTraits :: t_Object t_Object;
          typedef typename T_FACTORY :: t_Populator t_Populator;
          Operator::MakePopulator<t_Individual, t_Populator>::transform_unary
          ( 
            boost::bind( &random_two_valued<t_Object>, _1, _one, _two ),
            _function 
          );
          Print::xmg << Print::Xmg::comment << "Random" << Print::endl;
        }
     
    } // namespace AlloyLayers
  } // namespace GA
} // namespace LaDa
#endif
