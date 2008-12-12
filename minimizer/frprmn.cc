//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>

#include "frprmn.h"

namespace LaDa
{
  namespace Minimizer
  {
    namespace details
    {
      void* frpr_pointer_ = NULL;                           
    }

    inline bool Frpr :: Load_( const TiXmlElement &_element )
    {
      _element.Attribute( "itermax", &itermax );
      if ( not itermax ) itermax = 500;
      _element.Attribute( "tolerance", &tolerance );
      if ( not tolerance ) tolerance = types::tolerance;
      _element.Attribute( "linetolerance", &line_tolerance);
      if ( not line_tolerance ) line_tolerance = tolerance;
      _element.Attribute( "zeps", &zeps );
      if ( not zeps ) zeps = tolerance;
    
      return true;
    }


    inline const TiXmlElement* Frpr :: find_node( const TiXmlElement &_node )
    {
      const TiXmlElement *parent;
      std::string str;
    
      // This whole section tries to find a <Functional type="vff"> tag
      // in _element or its child
      str = _node.Value();
      if ( str.compare("Minimizer" ) != 0 )
        parent = _node.FirstChildElement("Minimizer");
      else
        parent = &_node;
    
      
      while (parent)
      {
        str = "";
        if ( parent->Attribute( "type" )  )
          str = parent->Attribute("type");
        if ( str.compare("frprmn" ) == 0 )
          break;
        parent = parent->NextSiblingElement("Minimizer");
      }
      __DOASSERT( not parent, 
                  "Could not find an <Minimizer type=\"frprmn\"> tag in input file.\n" )
      return parent;
    }

    bool Frpr :: Load( const TiXmlElement &_node )
    {
      const TiXmlElement* parent = find_node( _node );
      if( parent ) return Load_(*parent);
      return false;
    }

  }
} // namespace LaDa
