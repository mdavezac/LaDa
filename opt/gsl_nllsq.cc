//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "gsl_nllsq.h"
#include "gsl.h"


namespace Fitting
{
  bool NonLinearGsl :: Load( const TiXmlElement &_node )
  {
    const TiXmlElement *parent = &_node;
    for(; parent; parent = parent->NextSiblingElement("Fit") )
    {
      if( not parent->Attribute("type") ) continue;
      std::string name = parent->Attribute("type");
      if( name.compare( "Gsl" ) ) break;
    }
    if( not parent )
    {
      parent = _node.FirstChildElement("Fit" );
      if( not parent ) return false;
      return Load( *parent );
    }
    if( parent->Attribute( "tolerance" ) )
      parent->Attribute( "tolerance", &tolerance );
    types::t_int d( itermax );
    if( parent->Attribute( "itermax" ) )
      parent->Attribute( "itermax", &d );
    itermax = std::abs( d );
    return true;
  }

}
