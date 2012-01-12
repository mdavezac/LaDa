#include "LaDaConfig.h"

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL crystal_ARRAY_API
#include <numpy/arrayobject.h>

#include <iterator> 

#include <python/numpy_types.h>

#include "atom.h"
#include "structure_data.h"

namespace lp = LaDa::python;

using namespace LaDa;
using namespace LaDa::python;


#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

static PyMethodDef crystal_methods[] = {NULL};

PyMODINIT_FUNC initcppwrappers(void) 
{
  import_array(); // needed for NumPy 

  if (PyType_Ready(LaDa::crystal::atom_type()) < 0) return;
  if (PyType_Ready(LaDa::crystal::structure_type()) < 0) return;
  if (PyType_Ready(LaDa::crystal::structureiterator_type()) < 0) return;

  Py_INCREF(LaDa::crystal::atom_type());
  Py_INCREF(LaDa::crystal::structure_type());
  Py_INCREF(LaDa::crystal::structureiterator_type());

  char const doc[] =  "Wrapper around C++ atom/structure class and affiliates.";
  PyObject* module = Py_InitModule3("cppwrappers", crystal_methods, doc);

  PyModule_AddObject(module, "Atom", (PyObject *)LaDa::crystal::atom_type());
  PyModule_AddObject(module, "Structure", (PyObject *)LaDa::crystal::structure_type());
}
