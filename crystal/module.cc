#include "LaDaConfig.h"

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL lada_math_ARRAY_API
#include <numpy/arrayobject.h>

#include <iterator> 

#include <crystal/python/numpy_types.h>
#include <crystal/python/wrap_numpy.h>

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif
#define LADA_CRYSTAL_MODULE 0
#include "crystal.h"

#include "atom/pybase.cc"
#include "structure/pybase.cc"
#include "hart-forcade/pybase.cc"
#include "utilities.cc"
#include "map_sites.cc"
#include "equivalent_structures.cc"
#include "primitive.cc"
#include "space_group.cc"
#include "neighbors.cc"
#include "coordination_shells.cc"
#include "confsplit.cc"
#include "periodic_dnc.cc"

#include "methods.hpp"

PyMODINIT_FUNC initcppwrappers(void) 
{
  using namespace LaDa::crystal;
  static void *api_capsule[BOOST_PP_SLOT(1)];
  PyObject *c_api_object;

  char const doc[] =  "Wrapper around C++ atom/structure class and affiliates.";
  PyObject* module = Py_InitModule3("cppwrappers", methods_table, doc);
  if(not module) return;

  import_array(); // needed for NumPy 

  /* Initialize the C API pointer array */
# undef PYLADA_CRYSTALMODULE_H
# include "crystal.h"

  /* Create a Capsule containing the API pointer array's address */
  c_api_object = PyCapsule_New((void *)api_capsule, "lada.crystal.cppwrappers._C_API", NULL);

  if (c_api_object != NULL) PyModule_AddObject(module, "_C_API", c_api_object);

  if (PyType_Ready(atom_type()) < 0) return;
  if (PyType_Ready(structure_type()) < 0) return;
  if (PyType_Ready(structureiterator_type()) < 0) return;
  if (PyType_Ready(hftransform_type()) < 0) return;

  Py_INCREF(atom_type());
  Py_INCREF(structure_type());
  Py_INCREF(structureiterator_type());
  Py_INCREF(hftransform_type());

  PyModule_AddObject(module, "Atom", (PyObject *)atom_type());
  PyModule_AddObject(module, "Structure", (PyObject *)structure_type());
  PyModule_AddObject(module, "HFTransform", (PyObject *)hftransform_type());
}
