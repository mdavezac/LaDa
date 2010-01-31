//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/borrowed.hpp>

#include "vff.hpp"

BOOST_PYTHON_MODULE(vff)
{
  namespace bp = boost::python;
  bp::scope scope;
  scope.attr("__doc__") = "Valence Force-Field functional for the Zinc-Blende lattice.\n";
  bp::docstring_options doc_options(true, false);

  // loads lada.math first
  namespace bp = boost::python;
  bp::handle<> math( bp::borrowed(PyImport_ImportModule("lada.math")) );

  LaDa::Python::expose_vff();
  LaDa::Python::expose_layeredvff();
}
