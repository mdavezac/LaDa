#include "LaDaConfig.h"

#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>

#include <opt/types.h>
#include "redirect.hpp"


BOOST_PYTHON_MODULE(_opt)
{
  namespace bp = boost::python;
  bp::scope scope;
  bp::docstring_options doc_options(true, false);
  scope.attr("__doc__") = "imported into L{opt} namespace.";

  LaDa::python::expose_redirect();
}
