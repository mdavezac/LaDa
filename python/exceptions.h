#ifndef LADA_PYTHON_EXCEPTIONS_H
#define LADA_PYTHON_EXCEPTIONS_H
#include "LaDaConfig.h"

// #include <boost/python/object.hpp>
// #include <boost/python/import.hpp>
// #include <boost/python/exception_translator.hpp>
// #include <boost/exception/diagnostic_information.hpp>

#include <Python.h>
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
#   define LADA_PYERROR(EXCEPTION, MESSAGE)                                       \
      {                                                                           \
        PyObject* err_module = PyImport_ImportModule("lada.error");               \
        if(err_module)                                                            \
        {                                                                         \
          PyObject *err_result = PyObject_GetAttrString(err_module, #EXCEPTION);  \
          if(not err_result) Py_DECREF(err_module);                               \
          else                                                                    \
          {                                                                       \
            PyErr_SetString(err_result, MESSAGE);                                 \
            Py_DECREF(err_module);                                                \
            Py_DECREF(err_result);                                                \
          }                                                                       \
        }                                                                         \
      }
#   define LADA_PYTHROW(EXCEPTION, MESSAGE)                                       \
      {                                                                           \
        LADA_PYERROR(#EXCEPTION, MESSAGE);                                        \
        BOOST_THROW_EXCEPTION(error::EXCEPTION() << error::string(MESSAGE));      \
      }

    //! \def LADA_PYERROR(EXCEPTION, MESSAGE)
    //!      Raises a python exception with a formatted message, but no c++ exception.
    //!      For formatting, see PyErr_Format from the python C API.
    //!      EXCEPTION should be an unqualified declared in python/exceptions.h.
#   define LADA_PYERROR_FORMAT(EXCEPTION, MESSAGE, OTHER) \
      {                                                                           \
        PyObject* err_module = PyImport_ImportModule("lada.error");               \
        if(err_module)                                                            \
        {                                                                         \
          PyObject *err_result = PyObject_GetAttrString(err_module, #EXCEPTION);  \
          if(not err_result) Py_DECREF(err_module);                               \
          else                                                                    \
          {                                                                       \
            PyErr_Format(err_result, MESSAGE, OTHER);                             \
            Py_DECREF(err_module);                                                \
            Py_DECREF(err_result);                                                \
          }                                                                       \
        }                                                                         \
      }
  }
}
# endif 
