//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <iostream>
#include <algorithm>
#include <boost/lambda/lambda.hpp>

#include <opt/types.h>
#include <opt/convex_hull.h>
#include <opt/debug.h>

#include "misc.hpp"
#include "xml.hpp"
#include "structure.hpp"

namespace PythonLaDa
{
  namespace XML
  {
    template<> std::string nodename<Crystal::Structure>()  { return "Structure"; }
  }
  using namespace boost::python;
  void export_structure()
  {
  
    using namespace boost::python;
     class_< Crystal::Structure >( "Structure" )
       .def( init< Crystal::Structure >() )
       .def_readwrite( "cell",   &Crystal::Structure::cell )
       .def_readwrite( "atoms",  &Crystal::Structure::atoms )
       .def_readwrite( "energy", &Crystal::Structure::energy )
       .def_readwrite( "scale",  &Crystal::Structure::scale )
       .def_readwrite( "index", &Crystal::Structure::name )
       .def( "__str__",  &print<Crystal::Structure> ) 
       .def( "fromXML",  &XML::from<Crystal::Structure> )
       .def( "toXML",  &XML::to<Crystal::Structure> );
  }

}
