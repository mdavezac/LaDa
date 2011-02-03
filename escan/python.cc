#include "LaDaConfig.h"

#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/borrowed.hpp>

#include <numpy/ndarrayobject.h>

#include "escan.hpp"
#include "wavefunctions.hpp"

BOOST_PYTHON_MODULE(_escan)
{
  namespace bp = boost::python;
  bp::docstring_options doc_options(true, false);
  bp::scope scope;
  scope.attr("__doc__") = "This namespace is imported into lada.escan.\n";
  scope.attr("__docformat__") = "restructuredtext en";

  // loads lada.math first
  namespace bp = boost::python;
  bp::handle<> math( bp::borrowed(PyImport_ImportModule("lada.math")) );

  LaDa::python::expose_escan();
  LaDa::python::expose_wfns();
}
