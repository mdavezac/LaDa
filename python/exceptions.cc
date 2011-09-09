#include "LaDaConfig.h"

#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/borrowed.hpp>
#include <boost/python/def.hpp>
#include <boost/python/str.hpp>
#include <boost/python/exception_translator.hpp>

#include <root_exceptions.h>

namespace bp = boost::python;

BOOST_PYTHON_MODULE(exceptions)
{
  bp::scope scope;
  scope.attr("__docformat__") = "restructuredtext en";
  scope.attr("__doc__") = "This namespace is imported into lada.crystal.\n\n"
                          "It contains all C++ exception types.";
  bp::docstring_options doc_options(true, false);

  namespace bp = boost::python;
  char root_name[] = "lada.exceptions.root";
  char root_doc[] = "Root of all lada exceptions.";
  scope.attr("root") = bp::object( bp::handle<>( PyErr_NewExceptionWithDoc(root_name, root_doc, NULL, NULL) ) ); 
}
