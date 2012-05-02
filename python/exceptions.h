#ifndef LADA_PYTHON_EXCEPTIONS_H
#define LADA_PYTHON_EXCEPTIONS_H
#include "LaDaConfig.h"

#include <boost/python/object.hpp>
#include <boost/python/import.hpp>
#include <boost/python/exception_translator.hpp>
#include <boost/exception/diagnostic_information.hpp>

#include <root_exceptions.h>
#include <math/exceptions.h>

namespace LaDa
{
  namespace error
  {
    //! Attribute error thrown explicitely by lada.
    struct AttributeError: virtual root {};
    //! Key error thrown explicitely by lada.
    struct KeyError: virtual out_of_range, virtual root {};
    //! Value error thrown explicitely by lada.
    struct ValueError: virtual root {};
    //! Index error thrown explicitely by lada.
    struct IndexError: virtual root {};
    //! Argument error thrown explicitely by lada.
    struct TypeError: virtual root {};
    //! Not implemented error thrown explicitely by lada.
    struct NotImplementedError: virtual root {};
    //! Subclasses python's ImportError.
    struct ImportError: virtual root {};

    boost::python::object inline get_error(std::string const &_name)
    {
      namespace bp = boost::python;
      bp::object const module = bp::import("lada.error");
      bool const has_attr = PyObject_HasAttrString(module.ptr(), _name.c_str());
      return has_attr ? module.attr(_name.c_str()): module.attr("root");
    }
#   ifdef LADA_ERROR 
#     error LADA_ERROR already defined.
#   endif
#   define LADA_ERROR(TYPE)                                                  \
    inline void translate ## TYPE(::LaDa::error::TYPE const &_e)             \
    {                                                                        \
      if(PyErr_Occurred() != NULL) return;                                   \
      boost::python::object const error = get_error(#TYPE);                  \
      std::string message = boost::diagnostic_information(_e);               \
      message += "lada.error.";                                              \
      message += #TYPE;                                                      \
      message += " encountered.";                                            \
      PyErr_SetString(error.ptr(), message.c_str());                         \
    }
    LADA_ERROR(root);
    LADA_ERROR(input);
    LADA_ERROR(internal);
    LADA_ERROR(out_of_range);
    LADA_ERROR(infinite_loop);
    LADA_ERROR(KeyError);
    LADA_ERROR(ValueError);
    LADA_ERROR(IndexError);
    LADA_ERROR(TypeError);
    LADA_ERROR(NotImplementedError);
    LADA_ERROR(ImportError);
    LADA_ERROR(AttributeError);
    LADA_ERROR(math);
    LADA_ERROR(singular_matrix);
#   undef LADA_ERROR
    inline void bp_register()
    {
      namespace bp = boost::python;
      static bool is_first = true;
      if(not is_first) return;
      is_first = false;
#     define LADA_ERROR(TYPE) \
        bp::register_exception_translator< LaDa::error::TYPE>(&translate ## TYPE);
        LADA_ERROR(root);
        LADA_ERROR(input);
        LADA_ERROR(internal);
        LADA_ERROR(out_of_range);
        LADA_ERROR(infinite_loop);
        LADA_ERROR(KeyError);
        LADA_ERROR(ValueError);
        LADA_ERROR(IndexError);
        LADA_ERROR(TypeError);
        LADA_ERROR(NotImplementedError);
        LADA_ERROR(ImportError);
        LADA_ERROR(AttributeError);
        LADA_ERROR(math);
        LADA_ERROR(singular_matrix);
#     undef LADA_ERROR
    }

#   ifdef LADA_PYERROR
#     error LADA_PYERROR already  defined. 
#   endif
#   ifdef LADA_PYERROR_FORMAT
#     error LADA_PYERROR already  defined. 
#   endif
#   ifdef LADA_PYTHROW
#     error LADA_PYERROR already  defined. 
#   endif
    //! \def LADA_PYERROR(EXCEPTION, MESSAGE)
    //!      Raises a python exception with the interpreter, but no c++ exception.
    //!      EXCEPTION should be an unqualified declared in python/exceptions.h.
#   define LADA_PYERROR(EXCEPTION, MESSAGE) \
      PyErr_SetString(::LaDa::error::get_error(#EXCEPTION).ptr(), MESSAGE)
    //! \def LADA_PYERROR(EXCEPTION, MESSAGE)
    //!      Raises a python exception with a formatted message, but no c++ exception.
    //!      For formatting, see PyErr_Format from the python C API.
    //!      EXCEPTION should be an unqualified declared in python/exceptions.h.
#   define LADA_PYERROR_FORMAT(EXCEPTION, MESSAGE, OTHER) \
      PyErr_Format(::LaDa::error::get_error(#EXCEPTION).ptr(), MESSAGE, OTHER)
    //! \def LADA_PYTHROW(EXCEPTION, MESSAGE)
    //!      Raises a boost exception where EXCEPTION is stored as pyexcetp and MESSAGE as string.
    //!      EXCEPTION should be an unqualified declared in python/exceptions.h.
    //!      This macro makes it easy to catch all thrown python exceptions in a single statement.
#   define LADA_PYTHROW(EXCEPTION, MESSAGE)                                                   \
    {                                                                                         \
      PyObject * const exception = ::LaDa::error::get_error(#EXCEPTION).ptr();                \
      std::ostringstream sstr;                                                                \
      sstr << "In " << __FILE__ << "(" << __LINE__ << "): " << MESSAGE;                       \
      PyErr_SetString(exception, sstr.str().c_str());                                         \
      BOOST_THROW_EXCEPTION( ::LaDa::error::EXCEPTION()                                       \
                             << ::LaDa::error::string(sstr.str()) );                          \
    }
  }
}
# endif 
