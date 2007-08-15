#include <functional>
#include <algorithm>
#include <ext/algorithm>

#include "ce.h"
#include "functors.h"
#include "print_xmgrace.h"
#include "lamarck/atom.h"
#include "opt/va_minimizer.h"
#include "opt/bisection.h"
#include "concentration.h"

namespace CE
{

  bool Evaluator :: Load( const TiXmlElement &_node )
  {
    if ( not t_Base :: Load ( _node ) )
    {
      std::cerr << " Could  not Load SingleSite stuff "
                << std::endl << "Giving up" << std::endl;
      return false;
    }
    const TiXmlElement *functional_xml = _node.FirstChildElement("Functional");
    for(; functional_xml; functional_xml = functional_xml->NextSiblingElement("Functional") )
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
      return false;
    Ising_CE::Structure::lattice = VA_CE::Functional_Builder::lattice;
       
    add_equivalent_clusters();

    generate_functional(structure, &functional);

    return true;
  }


} // namespace CE


