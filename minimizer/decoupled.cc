//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "decoupled.h"

namespace LaDa
{
  namespace Minimizer
  {

    bool Any::Load( const TiXmlElement& _node )
    {
      const TiXmlElement *node = opt::find_node( _node, "Minimizer" );
      for(; node; node = node->NextSiblingElement( "Minimizer" ) )
      {
        DOASSERT( not node->Attribute("type"), "Found minimizer tag without type attribute.\n" )
        if( std::string( node->Attribute( "type" ) ) == "decoupled" ) break;
      }
      if( not node ) return false;
      __DOASSERT( not nod->Attribute( "mid" ), "Decoupled minimizer must have \"mid\" attribute.\n" )
      mid = boost::lexical_cast< size_t >( node->Attribute( "mid" ) );
      if( node->Attribute("tolerance") )
        tolerance = boost::lexical_cast< types::t_real >( node->Attribute( "tolerance" ) );
      if( node->Attribute("itermax") )
        itermax = boost::lexical_cast< types::t_int >( node->Attribute( "itermax" ) );

      __DOASSERT( not node->FirstSiblingElement( "Minimizer" ), 
                  "Decoupled minimizer must contain nested Minimizer tag.\n" )
      return any.Load( *node->FirstSiblingElement( "Minimizer" ) );
    }

  }
} 
