//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "any.h"

namespace LaDa
{
  namespace Minimizer
  {

    bool Any::Load( const TiXmlElement& _node )
    {
      const TiXmlElement *node = opt::find_node( _node, "Minimizer" );
      boost::shared_ptr< Gsl > gsl( new Gsl );
      boost::shared_ptr< Frpr > frpr( new Frpr );
      for(; node; node = node->NextSiblingElement( "Minimizer" ) )
      {
        if( not node->Attribute( "type" ) ) continue;
        if( gsl->Load( *node ) )
        { 
          std::cout << "Loaded Gsl minimizer.\n";
          gsl_.swap( gsl );
          return true; 
        }
        if( frpr->Load( *node ) )
        { 
          std::cout << "Loaded original VFF minimizer.\n";
          frpr_.swap( frpr ); 
          return true; 
        }
      }

      return false;
    }

  }
} 
