#include "LaDaConfig.h"

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL enumeration_ARRAY_API
#include <numpy/arrayobject.h>

#include <python/exceptions.h>
#include <math/fuzzy.h>
#include <python/numpy_types.h>
#include <python/object.h>
#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

namespace LaDa
{
  namespace enumeration
  {
    //! Checks whether a numpy array is an integer. 
    PyObject* is_integer(PyObject *_module, PyObject *_array)
    {
      if(not PyArray_Check(_array))
      {
        LADA_PYTHROW(TypeError, "input must be a numpy array.\n");
        return NULL;
      }
      if(PyArray_ISINTEGER(_array)) Py_RETURN_TRUE;
      python::Object iterator(PyArray_IterNew(_array));
      if(not iterator) return NULL;
      PyObject* i_iterator = iterator.borrowed();
      int const type = ((PyArrayObject*)_array)->descr->type_num;
#     ifdef LADA_IFTYPE
#       error LADA_IFTYPE already defined
#     endif
#     define LADA_IFTYPE(TYPENUM, TYPE)                                        \
        if(type == TYPENUM)                                                    \
        {                                                                      \
          while(PyArray_ITER_NOTDONE(i_iterator))                              \
          {                                                                    \
            TYPE const a = *((TYPE*) PyArray_ITER_DATA(i_iterator));           \
            if(not math::is_null(std::floor(a+1e-5) - a)) Py_RETURN_FALSE;     \
            PyArray_ITER_NEXT(i_iterator);                                     \
          }                                                                    \
        }
      LADA_IFTYPE(NPY_FLOAT, math::numpy::type<npy_float>::np_type)
      else LADA_IFTYPE(NPY_DOUBLE, math::numpy::type<npy_double>::np_type)
      else LADA_IFTYPE(NPY_LONGDOUBLE, math::numpy::type<npy_longdouble>::np_type)
#     undef LADA_WITH_DATA_TYPE
      Py_RETURN_TRUE;
    }
    //! Methods table for crystal module.
    static PyMethodDef methods_table[] = {
        {"is_integer",  is_integer, METH_O,
         "True if the numpy array is an integer.\n\n"
         "This method takes a single argument which *must* be a *numpy* array.\n"},
        {NULL, NULL, 0, NULL}        /* Sentinel */
    };
  }
}

PyMODINIT_FUNC initcppwrappers(void) 
{
  import_array(); // needed for NumPy 
  LaDa::error::bp_register();

  char const doc[] =  "Wrapper around C++ enumeration methods.";
  PyObject* module = Py_InitModule3("cppwrappers", LaDa::enumeration::methods_table, doc);
}
