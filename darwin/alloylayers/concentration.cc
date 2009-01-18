//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <print/stdout.h>
#include <print/manip.h>
#include <print/xmg.h>

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/string.hpp>

#include "concentration.h"

namespace LaDa
{
  namespace GA
  {
    namespace Keepers 
    {
      bool Concentration :: Load ( const TiXmlElement &_node )
      {
        if ( _node.Attribute("x") )
          x = 2e0 * boost::lexical_cast< types::t_real > ( _node.Attribute("x") ) - 1e0;
        if ( _node.Attribute("y") )
          y = 2e0 * boost::lexical_cast< types::t_real > ( _node.Attribute("y") ) - 1e0;
        return true;
      }
      bool Concentration :: Save( TiXmlElement &_node ) const
      {
        _node.SetAttribute("vbm", boost::lexical_cast<std::string>( x * 0.5e0 + 0.5e0 ) );
        _node.SetAttribute("cbm", boost::lexical_cast<std::string>( y * 0.5e0 + 0.5e0 ) );
        return true;
      }
    } // namespace Keepers
  } // namespace GA
}


