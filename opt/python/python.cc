//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>

#include <python/std_vector.hpp>
#include <opt/types.h>

#include "convexhull.impl.hpp"
#include "errortuple.hpp"

BOOST_PYTHON_MODULE(_opt)
{
  namespace bp = boost::python;
  bp::scope scope;
  scope.attr("__doc__") = "imported into L{opt} namespace.";
  bp::docstring_options doc_options(true, false);

  LaDa::Python::expose_errors();
  LaDa::Python::exposeConvexHull<boost::python::object>( "ConvexHull" );
  LaDa::Python::expose_vector<LaDa::types::t_real>( "cReals", "A stl container of real values." );
}
