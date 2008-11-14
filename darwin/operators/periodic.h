//
//  Version: $Id$
//
#ifndef _LADA_DARWIN_OPERATORS_PERIODIC_H_
#define _LADA_DARWIN_OPERATORS_PERIODIC_H_

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include<boost/lexical_cast.hpp>
#include<boost/bind.hpp>
#include<boost/ref.hpp>

#include<opt/types.h>
#include<print/xmg.h>

#include "../gencount.h"

namespace LaDa
{
  namespace GA
  {
    namespace Operator
    {
      template< class T_POPULATOR > 
        void periodic( T_POPULATOR &_populator,
                       boost::function<void(T_POPULATOR&)>& _function,
                       types::t_unsigned &_period,
                       const GenCount &_age )
        { if ( _age() % _period == 0 ) _function( _populator ); }
    } // namespace Operator

    namespace Factory
    {
      template< class T_FACTORY >
        void periodic( T_FACTORY &_factory,
                       boost::function<void( typename T_FACTORY::t_Populator& )>&
                         _function,
                       const std::string &_value,
                       const GenCount &_age )
        {
          typedef typename T_FACTORY::t_Populator t_Populator;
          types::t_unsigned 
            period( boost::lexical_cast< types::t_unsigned >(_value) );
          _function = boost::bind( &Operator::periodic<t_Populator>,
                                   _1, _function, period, boost::cref(_age) );
          Print::xmg << Print::Xmg::addtolast << " period = " << period << Print::endl;
        }
        
    } // end of Factory namespace.
  } // namespace GA
} // namespace LaDa

#endif

