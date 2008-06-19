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
    std::sort( structure.k_vecs.begin(), structure.k_vecs.end(),  Crystal::sort_kvec );
    
    const TiXmlElement *functional_xml = _node.FirstChildElement("Functional");
    for(; functional_xml;
          functional_xml = functional_xml->NextSiblingElement("Functional") )
    {
      if ( not functional_xml->Attribute("type") ) continue;
      std::string str =  functional_xml->Attribute("type");
      if ( str.compare("CE") == 0 ) break;
    }
    __DOASSERT( not functional_xml,
                "No <Functional type=\"CE\"> found.\n"
                "Giving up\n" )
    __DOASSERT( not t_Base::Load(*functional_xml),
                "Could not load CE functional from XML.\n" )
    Crystal::Structure::lattice = t_Base::lattice;
       
    add_equivalent_clusters();

    generate_functional(structure, &functional);
    
    return true;
  }

#ifdef _MPI
  void Darwin :: set_mpi( boost::mpi::communicator *_comm )
  {
    functional.get_functional1()->set_mpi( _comm );
    functional.get_functional2()->set_mpi( _comm );
  }
#endif

} // namespace CE


