//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
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

  template< class BULL >
  Crystal::Lattice& return_crystal_lattice( BULL & )
  {
    __DOASSERT( not Crystal::Structure::lattice,
                "Lattice pointer has not been set.\n" )
    return *Crystal::Structure::lattice; 
  }

  void expose_structure()
  {
    using namespace boost::python;
    class_< Crystal::Structure >( "Structure" )
      .def( init< Crystal::Structure >() )
      .def_readwrite( "cell",    &Crystal::Structure::cell )
      .def_readwrite( "atoms",   &Crystal::Structure::atoms )
      .def_readwrite( "k_vecs",   &Crystal::Structure::k_vecs )
      .def_readwrite( "energy",  &Crystal::Structure::energy )
      .def_readwrite( "scale",   &Crystal::Structure::scale )
      .def_readwrite( "index", &Crystal::Structure::name )
      .def( "__str__",  &print<Crystal::Structure> ) 
      .def( "fromXML",  &XML::from<Crystal::Structure> )
      .def( "toXML",  &XML::to<Crystal::Structure> )
      .def( "lattice", &return_crystal_lattice< Crystal::Structure >,
            return_value_policy<reference_existing_object>() );
    class_< Crystal::TStructure<std::string> >( "sStructure" )
      .def( init< Crystal::TStructure<std::string> >() )
      .def_readwrite( "cell",    &Crystal::TStructure<std::string>::cell )
      .def_readwrite( "atoms",   &Crystal::TStructure<std::string>::atoms )
      .def_readwrite( "energy",  &Crystal::TStructure<std::string>::energy )
      .def_readwrite( "scale",   &Crystal::TStructure<std::string>::scale )
      .def_readwrite( "index", &Crystal::TStructure<std::string>::name )
      .def( "__str__",  &print< Crystal::TStructure<std::string> > ) 
      .def( "fromXML",  &XML::from< Crystal::TStructure<std::string> > )
      .def( "toXML",  &XML::to< Crystal::TStructure<std::string> > )
      .def( "lattice", &return_crystal_lattice< Crystal::TStructure<std::string> >,
            return_value_policy<reference_existing_object>() );
  }

}
