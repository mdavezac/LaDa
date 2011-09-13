#include "LaDaConfig.h"

#include <root_exceptions.h>

#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>

#include "exceptions.h"

namespace bp = boost::python;

BOOST_PYTHON_MODULE(error)
{
  bp::scope scope;
  scope.attr("__docformat__") = "restructuredtext en";
  scope.attr("__doc__") = "This namespace is imported into lada.crystal.\n\n"
                          "It contains all C++ exception types.";
  bp::docstring_options doc_options(true, false);

  LADA_REGISTER_PYEXCEPT( ::LaDa::error::root, "lada.error.root",
                          "Root of all exceptions explicitely thrown by lada.", scope );
  LADA_REGISTER_PYEXCEPT_WITH_BASE( ::LaDa::error::input, "lada.error.input",
                          "Root of all input exceptions explicitely thrown by lada.", scope,
                          bp::object(bp::handle<>(PyExc_ValueError)) );
  LADA_REGISTER_PYEXCEPT( ::LaDa::error::internal, "lada.error.internal",
                          "Root of all internal (bugs) exceptions explicitely thrown by lada.", scope);
  LADA_REGISTER_PYEXCEPT_WITH_BASE( ::LaDa::error::out_of_range, "lada.error.out_of_range",
                          "Root of all out_of_range exceptions explicitely thrown by lada.", scope,
                          bp::object(bp::handle<>(PyExc_IndexError)));
  LADA_REGISTER_PYEXCEPT( ::LaDa::error::infinite_loop, "lada.error.infinite_loop",
                          "Root of all exceptions thrown by lada when it encounters an infinite loop.", scope );
}
