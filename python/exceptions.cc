#include "LaDaConfig.h"

#include <root_exceptions.h>

#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/borrowed.hpp>

#include "exceptions.h"

namespace bp = boost::python;

BOOST_PYTHON_MODULE(error)
{
  using namespace ::LaDa::error;
  using namespace ::LaDa::python;
  bp::scope scope;
  scope.attr("__docformat__") = "restructuredtext en";
  scope.attr("__doc__") = "This namespace is imported into lada.crystal.\n\n"
                          "It contains all C++ exception types.";
  bp::docstring_options doc_options(true, false);

  LADA_REGISTER_PYEXCEPT( root, "lada.error.root",
                          "Root of all exceptions explicitely thrown by lada.", scope );
  LADA_REGISTER_PYEXCEPT_WITH_BASE( input, "lada.error.input",
                          "Root of all input exceptions explicitely thrown by lada.", scope,
                          bp::make_tuple(bp::object(bp::handle<>(PyExc_ValueError))) );
  LADA_REGISTER_PYEXCEPT( internal, "lada.error.internal",
                          "Root of all internal (bugs) exceptions explicitely thrown by lada.", scope);
  LADA_REGISTER_PYEXCEPT_WITH_BASE( ::LaDa::error::out_of_range, "lada.error.out_of_range",
                          "Root of all out_of_range exceptions explicitely thrown by lada.", scope,
                          bp::make_tuple(bp::object(bp::handle<>(PyExc_IndexError))));
  LADA_REGISTER_PYEXCEPT( infinite_loop, "lada.error.infinite_loop",
                          "Root of all exceptions thrown by lada when it encounters an infinite loop.", scope );

  LADA_REGISTER_PYEXCEPT_WITH_BASE( ValueError, "lada.error.ValueError",
                          "Value errors explicitely thrown by lada.", scope,
                          bp::make_tuple( PyException<input>::exception(), 
                                          bp::object(bp::borrowed<>(PyExc_ValueError))));
  LADA_REGISTER_PYEXCEPT_WITH_BASE( KeyError, "lada.error.KeyError",
                          "Key errors explicitely thrown by lada.", scope,
                          bp::make_tuple( PyException<input>::exception(), 
                                          PyException<out_of_range>::exception(),
                                          bp::object(bp::borrowed<>(PyExc_KeyError))));
  LADA_REGISTER_PYEXCEPT_WITH_BASE( AttributeError, "lada.error.AttributeError",
                          "Attribute errors explicitely thrown by lada.", scope,
                          bp::make_tuple( PyException<input>::exception(),
                                          bp::object(bp::borrowed<>(PyExc_AttributeError))));
  LADA_REGISTER_PYEXCEPT_WITH_BASE( IndexError, "lada.error.IndexError",
                          "Index errors explicitely thrown by lada.", scope,
                          bp::make_tuple( PyException<input>::exception(),
                                          PyException<out_of_range>::exception(),
                                          bp::object(bp::borrowed<>(PyExc_IndexError))));
  LADA_REGISTER_PYEXCEPT_WITH_BASE( TypeError, "lada.error.TypeError",
                          "Argument errors explicitely thrown by lada.", scope,
                          bp::make_tuple( PyException<input>::exception(), 
                                          bp::object(bp::borrowed<>(PyExc_TypeError))) );
  LADA_REGISTER_PYEXCEPT_WITH_BASE( NotImplementedError, "lada.error.NotImplementedError",
                          "Stuff that has not yet been implemented.", scope,
                          bp::make_tuple( PyException<internal>::exception(),
                                          bp::object(bp::borrowed<>(PyExc_NotImplementedError))));
  LADA_REGISTER_PYEXCEPT_WITH_BASE( InternalError, "lada.error.InternalError",
                          "Internal error.", scope,
                          bp::make_tuple( PyException<internal>::exception(),
                                          bp::object(bp::borrowed<>(PyExc_RuntimeError))));
}
