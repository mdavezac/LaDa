#include "LaDaConfig.h"

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL enumeration_ARRAY_API
#include <numpy/arrayobject.h>

#include <python/exceptions.h>
#include <math/fuzzy.h>
#include <crystal/python/numpy_types.h>
#include <crystal/python/object.h>
#include "ndimiterator.h"
#include "fciterator.h"
#include "manipulations.h"
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
        LADA_PYERROR(TypeError, "input must be a numpy array.\n");
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
    PyObject* lexcompare(PyObject *_module, PyObject *_args)
    {
#     ifdef LADA_DEBUG
      if(_args == NULL)
      {
        LADA_PYERROR(TypeError, "_lexcompare expects two arguments.");
        return NULL;
      }
      if(PyTuple_Size(_args) != 2)
      {
        LADA_PYERROR(TypeError, "_lexcompare expects two arguments.");
        return NULL;
      }
#     endif
      PyArrayObject *first = (PyArrayObject*) PyTuple_GET_ITEM(_args, 0);
      PyArrayObject *second = (PyArrayObject*) PyTuple_GET_ITEM(_args, 1);
#     ifdef LADA_DEBUG
      if(not PyArray_Check(first))
      {
        LADA_PYERROR(TypeError, "First argument to _lexcompare is not a numpy array.");
        return NULL;
      }
      if(not PyArray_Check(second))
      {
        LADA_PYERROR(TypeError, "Second argument to _lexcompare is not a numpy array.");
        return NULL;
      }
      if(PyArray_NDIM(first) != 1 or PyArray_NDIM(second) != 1)
      {
        LADA_PYERROR(TypeError, "_lexcompare arguments should be 1d array.");
        return NULL;
      }
      if(PyArray_DIM(first, 0) != PyArray_DIM(second, 0))
      {
        LADA_PYERROR(TypeError, "_lexcompare arguments should have the same size.");
        return NULL;
      }
      if( first->descr->type_num != math::numpy::type<t_ndim>::value
          or second->descr->type_num != math::numpy::type<t_ndim>::value )
      {
        LADA_PYERROR(TypeError, "Wrong kind for _lexcompare arguments.");
        return NULL;
      }
#     endif

      npy_intp const stridea = PyArray_STRIDE(first, 0);
      npy_intp const strideb = PyArray_STRIDE(second, 0);
      npy_intp const n = PyArray_DIM(first, 0);
      char * dataa = (char*) PyArray_DATA(first);
      char * datab = (char*) PyArray_DATA(second);
      for(npy_intp i(0); i < n; ++i, dataa += stridea, datab += strideb)
      {
        t_ndim const a = *((t_ndim*) dataa);
        t_ndim const b = *((t_ndim*) datab);
        if(a > b) { return PyInt_FromLong(1); }
        else if(a < b) { return PyInt_FromLong(-1); }
      }
      return PyInt_FromLong(0); 
    }
    //! Methods table for crystal module.
    static PyMethodDef methods_table[] = {
        {"is_integer",  is_integer, METH_O,
         "True if the numpy array is an integer.\n\n"
         "This method takes a single argument which *must* be a *numpy* array.\n"},
        {"_lexcompare",  lexcompare, METH_VARARGS,
         "Lexicographic compare of two numpy arrays.\n\n"
         "Read the code for this function. If you do not understand it, do not\n"
         "use it.\n\n"
         ":returns:\n\n"
         "   - a > b: 1\n"
         "   - a == b: 0\n"
         "   - a < b: -1\n" },
        {NULL, NULL, 0, NULL}        /* Sentinel */
    };
  }
}

PyMODINIT_FUNC initcppwrappers(void) 
{
  import_array(); // needed for NumPy 

  if (PyType_Ready(LaDa::enumeration::ndimiterator_type()) < 0) return;
  Py_INCREF(LaDa::enumeration::ndimiterator_type());
  if (PyType_Ready(LaDa::enumeration::manipulations_type()) < 0) return;
  Py_INCREF(LaDa::enumeration::manipulations_type());
  if (PyType_Ready(LaDa::enumeration::fciterator_type()) < 0) return;
  Py_INCREF(LaDa::enumeration::fciterator_type());

  char const doc[] =  "Wrapper around C++ enumeration methods.";
  PyObject* module = Py_InitModule3("cppwrappers", LaDa::enumeration::methods_table, doc);
  PyModule_AddObject(module, "NDimIterator", (PyObject *)LaDa::enumeration::ndimiterator_type());
  PyModule_AddObject(module, "Manipulations", (PyObject *)LaDa::enumeration::manipulations_type());
  PyModule_AddObject(module, "FCIterator", (PyObject *)LaDa::enumeration::fciterator_type());
}
