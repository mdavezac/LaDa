#ifndef PYLADA_PYTHON_EXCEPTIONS_H
#define PYLADA_PYTHON_EXCEPTIONS_H
#include "PyladaConfig.h"

#include <Python.h>
#include <root_exceptions.h>

namespace Pylada
{
  namespace error
  {
    //! Attribute error thrown explicitely by pylada.
    struct AttributeError: virtual root {};
    //! Key error thrown explicitely by pylada.
    struct KeyError: virtual out_of_range, virtual root {};
    //! Value error thrown explicitely by pylada.
    struct ValueError: virtual root {};
    //! Index error thrown explicitely by pylada.
    struct IndexError: virtual root {};
    //! Argument error thrown explicitely by pylada.
    struct TypeError: virtual root {};
    //! Not implemented error thrown explicitely by pylada.
    struct NotImplementedError: virtual root {};
    //! Subclasses python's ImportError.
    struct ImportError: virtual root {};

#   ifdef PYLADA_PYERROR
#     error PYLADA_PYERROR already  defined. 
#   endif
#   ifdef PYLADA_PYERROR_FORMAT
#     error PYLADA_PYERROR already  defined. 
#   endif
#   ifdef PYLADA_PYTHROW
#     error PYLADA_PYERROR already  defined. 
#   endif
    //! \def PYLADA_PYERROR(EXCEPTION, MESSAGE)
    //!      Raises a python exception with the interpreter, but no c++ exception.
    //!      EXCEPTION should be an unqualified declared in python/exceptions.h.
#   define PYLADA_PYERROR(EXCEPTION, MESSAGE)                                       \
      {                                                                           \
        PyObject* err_module = PyImport_ImportModule("pylada.error");               \
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
#   define PYLADA_PYTHROW(EXCEPTION, MESSAGE)                                       \
      {                                                                           \
        PYLADA_PYERROR(EXCEPTION, MESSAGE);                                         \
        BOOST_THROW_EXCEPTION(error::EXCEPTION() << error::string(MESSAGE));      \
      }

    //! \def PYLADA_PYERROR(EXCEPTION, MESSAGE)
    //!      Raises a python exception with a formatted message, but no c++ exception.
    //!      For formatting, see PyErr_Format from the python C API.
    //!      EXCEPTION should be an unqualified declared in python/exceptions.h.
#   define PYLADA_PYERROR_FORMAT(EXCEPTION, MESSAGE, OTHER) \
      {                                                                           \
        PyObject* err_module = PyImport_ImportModule("pylada.error");               \
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
