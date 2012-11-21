#include "LaDaConfig.h"

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL lada_ARRAY_API
#include <numpy/arrayobject.h>

#include <errors/exceptions.h>
// #include <crystal/python/numpy_types.h>
// #include <crystal/python/wrap_numpy.h>

#define LADA_MATH_MODULE

#include "fuzzy.h"
#include "gruber.h"
#include "smith_normal_form.h"
#include "methods.hpp"
#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

PyMODINIT_FUNC initmath(void) 
{
  char const doc[] =  "A ratatouille of useful mathematics.";
  PyObject* module = Py_InitModule3("math", LaDa::math::methods_table, doc);
  if(not module) return;
  import_array(); // needed for NumPy 
}
