#include "LaDaConfig.h"

#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/borrowed.hpp>

#include "clj.hpp"
#include "functional.hpp"

BOOST_PYTHON_MODULE(pcm)
{
  namespace bp = boost::python;
  bp::scope scope;
  scope.attr("__doc__") = "Point Charge Model functional.";
  bp::docstring_options doc_options(true, false);

  // loads lada.math first
  namespace bp = boost::python;
  bp::handle<> math( bp::borrowed(PyImport_ImportModule("lada.math")) );

  LaDa::python::expose_clj();
  LaDa::python::expose_functional();
}
