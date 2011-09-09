#include "LaDaConfig.h"

#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/borrowed.hpp>
#include <boost/python/def.hpp>
#include <boost/python/str.hpp>
#include <boost/python/exception_translator.hpp>
#include <boost/exception/get_error_info.hpp>

#include <root_exceptions.h>
#include "exceptions.h"

namespace bp = boost::python;

namespace LaDa
{
  namespace python
  {
    boost::python::object root_exception()
    {
      char root_name[] = "lada.exceptions.root";
      char root_doc[] = "Root of all lada exceptions.";
      static boost::python::object object( bp::handle<>(
            PyErr_NewExceptionWithDoc(root_name, root_doc, NULL, NULL) ) ); 
      return object;
    }

    void root_translate(error::root const &_e)
    {
      std::string message = _e.what() ;
      if(std::string const *str = boost::get_error_info<error::string>(_e))
        message += "\n" + (*str);
      message += "\n";
      PyErr_SetString(root_exception().ptr(), message.c_str());
    }
  }
}

void nothrow () {};
template<class T> void dothrow_nomessage() { BOOST_THROW_EXCEPTION(T()); }
template<class T> void dothrow_message()
{ BOOST_THROW_EXCEPTION(T() << LaDa::error::string("This is a message.")); }


BOOST_PYTHON_MODULE(exceptions)
{
  bp::scope scope;
  scope.attr("__docformat__") = "restructuredtext en";
  scope.attr("__doc__") = "This namespace is imported into lada.crystal.\n\n"
                          "It contains all C++ exception types.";
  bp::docstring_options doc_options(true, false);

  scope.attr("root") = LaDa::python::root_exception(); // initialize root exceptions
  bp::register_exception_translator<LaDa::error::root>(&LaDa::python::root_translate);
}
