//
//  Version: $Id$
//
#include <functional>
#include <algorithm>
#include <ext/algorithm>

#include "opt/va_minimizer.h"

#include "ce.h"
#include "functors.h"
#include "lamarck/atom.h"
#include "concentration.h"

namespace CE
{
  bool Evaluator :: Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type )
  {
    double d;
    if ( not _node.Attribute("enthalpy", &d ) ) goto errorout;
    _indiv.quantities() = (types::t_real) d;

    return t_Base::Load( _indiv, _node, _type );
errorout:
    std::cerr << "Could not Load CE::Object" << std::endl;
    return false;
  }
  bool Evaluator :: Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const
  { 
    double d = (double) _indiv.const_quantities();
    _node.SetDoubleAttribute("enthalpy", d );
    return t_Base::Save( _indiv, _node, GA::LOADSAVE_SHORT );
  }

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


