#include "LaDaConfig.h"

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL lada_ewald_ARRAY_API
#include <numpy/arrayobject.h>

#include <algorithm>

#include <errors/exceptions.h>
#include <python/python.h>
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
  if(not module) return;
  import_array();
  if(not LaDa::python::import()) return;
}
