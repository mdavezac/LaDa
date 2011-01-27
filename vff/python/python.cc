#include "LaDaConfig.h"

#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/borrowed.hpp>

#include "vff.hpp"
#include "data.hpp"

BOOST_PYTHON_MODULE(_vff)
{
  namespace bp = boost::python;
  bp::scope scope;
  scope.attr("__doc__") = "Valence Force-Field functional for the Zinc-Blende lattice.\n";
  bp::docstring_options doc_options(true, false);

  // loads lada.math first
  namespace bp = boost::python;
  bp::handle<> math( bp::borrowed(PyImport_ImportModule("lada.math")) );

  LaDa::python::expose_vff();
  LaDa::python::expose_layeredvff();
  LaDa::python::expose_data();
}
