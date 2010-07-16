#include "LaDaConfig.h"

#include <boost/lexical_cast.hpp>
#include "gsl_nllsq.h"
#include "gsl.h"


namespace LaDa
{
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
      itermax = 40;
      if( parent->Attribute( "itermax" ) )
        itermax = boost::lexical_cast<types::t_unsigned>( parent->Attribute("itermax") ); 
      return true;
    }

  }
} // namespace LaDa
