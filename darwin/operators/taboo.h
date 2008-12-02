//
//  Version: $Id: taboos.h 856 2008-11-14 17:00:23Z davezac $
//
#ifndef _LADA_GA_OPERATOR_TABOO_H_
#define _LADA_GA_OPERATOR_TABOO_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../taboos/container.h"

namespace LaDa
{
  namespace GA
  {
    namespace Operator
    {
      template< class T_INDIVIDUAL, class T_POPULATOR >
        void taboo( T_POPULATOR &_populator, 
                    Taboo::Container<T_INDIVIDUAL> &_taboo,
                    types::t_unsigned _max,
                    const boost::function<void(T_POPULATOR& )>& _inner,
                    const boost::function<void(T_POPULATOR& )>& _random )
        {
          types::t_unsigned  i = 0;
          do
          {
            ++i;
            _inner( _populator );
            if ( not _taboo( *_populator ) ) return;
            *_populator = _populator.select();
          } while ( i < _max );
      
          do 
          {
            ++i;
            _random( _populator );
            if ( not _taboo( *_populator ) ) return;
          } while ( i < UINT_MAX );
      
          std::ostringstream sstr;
          sstr << __LINE__ << ", line: " << __LINE__ << "\n"
               << "Could not create a non-taboo individual\n";
          throw std::runtime_error( sstr.str() );
        }

      template< class T_FACTORY >
        void taboo_factory( T_FACTORY &_factory,
                            boost::function<void( typename T_FACTORY::t_Populator& )>&
                              _function,
                            const TiXmlElement &_node,
                            Taboo::Container<typename T_FACTORY::t_Individual> &_taboo,
                            const std::string &_random )
        {
          typedef typename T_FACTORY::t_Individual t_Individual;
          typedef typename T_FACTORY::t_Populator t_Populator;
          typedef boost::function<void(t_Populator&)> t_Function;
       
          types::t_unsigned max(100);
          if( _node.Attribute( "max" ) )
            max = boost::lexical_cast< types::t_unsigned >( _node.Attribute("max") );
          Print::xmg << Print::Xmg::comment << "TabooOperator( random fallback: " 
                     << _random << " )" << Print::endl << Print::Xmg::indent;
          t_Function random;
          // creates random operator.
          _factory( _random, random, _node );
          Print::xmg << Print::Xmg::removelast;
          //! creates inner operator.
          _factory( _function, _node );
       
          // binds all.
          _function = boost::bind
                      (
                        &taboo<t_Individual, t_Populator>,
                        _1, boost::ref(_taboo), max, _function, random
                      );
          Print::xmg << Print::Xmg::unindent;
        }
    } // namespace Operator
  } // namespace GA
} // namespace LaDa 


#endif
