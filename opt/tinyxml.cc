//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "tinyxml.h"

namespace opt
{
  //! \brief Returns the node \a _name.
  //! \details Looks first to \a _element, then its childrent, then its
  //!          next siblings.
  //! \todo Look amongst all siblings.
  const TiXmlElement* find_node( const TiXmlElement &_element,
                                 const std::string& _name )
  {
    const TiXmlElement *parent;
  
    // Find first XML "Structure" node (may be _element from start, or a child of _element)
    std::string str = _element.Value();
    if ( str.compare(_name) == 0 ) return &_element;
    parent =  _element.FirstChildElement(_name);
    if( parent ) return parent;

    return _element.NextSiblingElement(_name);
  }


  //! \brief Returns the node \<Functional type=\a _name\>.
  //! \details Looks first to \a _element, then its childrent, then its
  //!          next siblings.
  const TiXmlElement* find_functional_node ( const TiXmlElement &_element,
                                             const std::string &_name )
  {
    const TiXmlElement *parent;
    std::string str;

    // This whole section tries to find a <Functional type="vff"> tag
    // in _element or its child
    str = _element.Value();
    if ( str.compare("Functional" ) != 0 )
      parent = _element.FirstChildElement("Functional");
    else parent = &_element;
    
    while (parent)
    {
      str = "";
      if ( parent->Attribute( "type" )  )
        str = parent->Attribute("type");
      if ( str.compare(_name) == 0 )
        break;
      parent = parent->NextSiblingElement("Functional");
    }
    if ( parent ) return parent;
    
    std::cerr << "Could not find an <Functional type=\"" << _name << "\"> tag in input file" 
              << std::endl;
    return NULL;
  }  // Functional :: find_node
}

