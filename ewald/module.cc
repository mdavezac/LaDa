#include "LaDaConfig.h"

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL lada_math_ARRAY_API
#include <numpy/arrayobject.h>

#include <algorithm>

#include <python/exceptions.h>
#include <python/numpy_types.h>
#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

#include "ewald.h"

namespace LaDa
{
  namespace pcm
  {
    //! Methods table for crystal module.
    static PyMethodDef methods_table[] = {
        {"ewald", (PyCFunction)ewald, METH_KEYWORDS,
         "Performs an ewald summation." },
        {NULL, NULL, 0, NULL}        /* Sentinel */
    };
  }
}

PyMODINIT_FUNC initcppwrappers(void) 
{
  import_array(); // needed for NumPy 

  char const doc[] =  "Wrapper around C++/fortan point-ion models methods.";
  PyObject* module = Py_InitModule3("cppwrappers", LaDa::pcm::methods_table, doc);
}
