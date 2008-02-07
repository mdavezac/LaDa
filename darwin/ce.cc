//
//  Version: $Id$
//
#include <functional>
#include <algorithm>

#include <lamarck/atom.h>
#include <lamarck/structure.h>

#include "ce.h"
#include "functors.h"
#include "concentration.h"

namespace CE
{
  Darwin::~Darwin() 
  {
    if ( functional.get_functional1() ) 
      delete functional.get_functional1();
    if ( functional.get_functional2() ) 
      delete functional.get_functional2();
  }


  bool Darwin :: Load( const TiXmlElement &_node )
  {
    std::sort( structure.k_vecs.begin(), structure.k_vecs.end(),  Ising_CE::sort_kvec );
    
    const TiXmlElement *functional_xml = _node.FirstChildElement("Functional");
    for(; functional_xml;
          functional_xml = functional_xml->NextSiblingElement("Functional") )
    {
      std::string str = ""; 
      if ( functional_xml->Attribute("type") )
        str = functional_xml->Attribute("type");
      if ( str.compare("CE") == 0 )
        break;
    }
    if ( not functional_xml )
    {
      std::cerr << "No <Functional type=\"CE\"> found "
                << std::endl << "Giving up" << std::endl;
      return false;
    }

    if ( not VA_CE::Functional_Builder::Load(*functional_xml) )
    {
      std::cerr << "Could not load CE functional from XML\n";
      return false;
    }
    Ising_CE::Structure::lattice = VA_CE::Functional_Builder::lattice;
       
    add_equivalent_clusters();

    generate_functional(structure, &functional);
    
    return true;
  }

} // namespace CE


