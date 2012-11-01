#include "LaDaConfig.h"

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL lada_ARRAY_API
#include <numpy/arrayobject.h>

#include <python/numpy_types.h>
#include <python/wrap_numpy.h>

#include "fuzzy.h"
#include "methods.hpp"
#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

PyMODINIT_FUNC initmath(void) 
{
  LaDa::error::bp_register();

  char const doc[] =  "Wrapper around C++ atom/structure class and affiliates.";
  PyObject* module = Py_InitModule3("math", LaDa::math::methods_table, doc);
  import_array(); // needed for NumPy 
}
