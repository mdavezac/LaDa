//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <boost/python.hpp>

#include <opt/types.h>
#include <opt/debug.h>

#include "misc.hpp"
#include "xml.hpp"
#include "structure.hpp"

namespace PythonLaDa
{
  namespace XML
  {
    template<> std::string nodename<Crystal::Structure>()  { return "Structure"; }
    template<> std::string nodename< Crystal::TStructure<std::string> >() 
      { return nodename<Crystal::Structure>(); }
  }
  void expose_structure()
  {
  
    using namespace boost::python;
    class_< Crystal::Structure >( "Structure" )
      .def( init< Crystal::Structure >() )
      .def_readwrite( "cell",    &Crystal::Structure::cell )
      .def_readwrite( "lattice", Crystal::Structure::lattice )
      .def_readwrite( "atoms",   &Crystal::Structure::atoms )
      .def_readwrite( "energy",  &Crystal::Structure::energy )
      .def_readwrite( "scale",   &Crystal::Structure::scale )
      .def_readwrite( "index", &Crystal::Structure::name )
      .def( "__str__",  &print<Crystal::Structure> ) 
      .def( "fromXML",  &XML::from<Crystal::Structure> )
      .def( "toXML",  &XML::to<Crystal::Structure> );
    class_< Crystal::TStructure<std::string> >( "sStructure" )
      .def( init< Crystal::TStructure<std::string> >() )
      .def_readwrite( "cell",    &Crystal::TStructure<std::string>::cell )
      .def_readwrite( "atoms",   &Crystal::TStructure<std::string>::atoms )
      .def_readwrite( "energy",  &Crystal::TStructure<std::string>::energy )
      .def_readwrite( "scale",   &Crystal::TStructure<std::string>::scale )
      .def_readwrite( "index", &Crystal::TStructure<std::string>::name )
      .def( "__str__",  &print< Crystal::TStructure<std::string> > ) 
      .def( "fromXML",  &XML::from< Crystal::TStructure<std::string> > )
      .def( "toXML",  &XML::to< Crystal::TStructure<std::string> > );
  }

}
