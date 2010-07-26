#include "LaDaConfig.h"

#include <boost/lexical_cast.hpp>

#include "decoupled.h"

namespace LaDa
{
  namespace Minimizer
  {
    LADA_REGISTER_MINIMIZER_VARIANT_SOURCE( Decoupled, "decoupled" )

    bool Decoupled::Load( const TiXmlElement& _node )
    {
      const TiXmlElement *node = opt::find_node( _node, "Minimizer" );
      for(; node; node = node->NextSiblingElement( "Minimizer" ) )
      {
        LADA_DO_NASSERT( not node->Attribute("type"), "Found minimizer tag without type attribute.\n" )
        if( std::string( node->Attribute( "type" ) ) == "decoupled" ) break;
      }
      if( not node ) return false;
      LADA_DO_NASSERT( not node->Attribute( "mid" ),
                  "Decoupled minimizer must have \"mid\" attribute.\n" )
      mid = boost::lexical_cast< size_t >( node->Attribute( "mid" ) );
      LADA_DO_NASSERT( mid == 0, "Attribute \"mid\" cannot be zero.\n" )
      if( node->Attribute("tolerance") )
        tolerance = boost::lexical_cast< types::t_real >( node->Attribute( "tolerance" ) );
      if( node->Attribute("itermax") )
        itermax = boost::lexical_cast< types::t_int >( node->Attribute( "itermax" ) );

      LADA_DO_NASSERT( not node->FirstChildElement( "Minimizer" ), 
                  "Decoupled minimizer must contain nested Minimizer tag.\n" )
      return minimizer.Load( *node->FirstChildElement( "Minimizer" ) );
    }

  }
} 
