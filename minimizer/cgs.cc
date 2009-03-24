//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/lexical_cast.hpp>
#include <opt/tinyxml.h>

#include "cgs.h"

namespace LaDa
{
  namespace Fitting
  {

    bool Cgs :: load( const TiXmlElement &_node )
    {
      const TiXmlElement *parent( opt::find_node( _node, "Minimizer", "type", "cgs" ) );
      if( not parent ) return false;
      if( parent->Attribute("tolerance") )
        tolerance = boost::lexical_cast< types::t_real >( parent->Attribute("tolerance") );
      if( parent->Attribute("itermax") )
        itermax = boost::lexical_cast< types::t_unsigned >( parent->Attribute("itermax") );
      if( _node.Attribute("verbose") ) 
      {
        const std::string value( _node.Attribute("verbose") );
        if( value == "true" or value == "TRUE" or value == "T" or value == "t" ) 
          verbose = true;
        else if( value == "false" or value == "FALSE" or value == "F" or value == "f" ) 
          verbose = false;
        else verbose = boost::lexical_cast<bool>( _node.Attribute("verbose") );
      }
      return true;
    }
  }
} // namespace LaDa
