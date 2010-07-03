#include "LaDaConfig.h"

#include <opt/types.h>

#include "frprmn.h"

namespace LaDa
{
  namespace Minimizer
  {
    LADA_REGISTER_MINIMIZER_VARIANT_SOURCE( Frpr, "original VFF" )
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
