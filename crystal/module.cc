#include "LaDaConfig.h"

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL lada_crystal_ARRAY_API
#include <numpy/arrayobject.h>

#include <iterator> 

#include <python/numpy_types.h>
#include <python/wrap_numpy.h>

#include "supercell.h"
#include "map_sites.h"
#include "equivalent_structures.h"
#include "primitive.h"
#include "space_group.h"
#include "neighbors.h"
#include "coordination_shells.h"
#include "confsplit.h"
#include "periodic_dnc.h"
#include "methods.hpp"
#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

PyMODINIT_FUNC initcppwrappers(void) 
{
  LaDa::error::bp_register();

  if (PyType_Ready(LaDa::crystal::atom_type()) < 0) return;
  if (PyType_Ready(LaDa::crystal::structure_type()) < 0) return;
  if (PyType_Ready(LaDa::crystal::structureiterator_type()) < 0) return;
  if (PyType_Ready(LaDa::crystal::hftransform_type()) < 0) return;

  Py_INCREF(LaDa::crystal::atom_type());
  Py_INCREF(LaDa::crystal::structure_type());
  Py_INCREF(LaDa::crystal::structureiterator_type());
  Py_INCREF(LaDa::crystal::hftransform_type());

  char const doc[] =  "Wrapper around C++ atom/structure class and affiliates.";
  PyObject* module = Py_InitModule3("cppwrappers", LaDa::crystal::methods_table, doc);
  import_array(); // needed for NumPy 

  PyModule_AddObject(module, "Atom", (PyObject *)LaDa::crystal::atom_type());
  PyModule_AddObject(module, "Structure", (PyObject *)LaDa::crystal::structure_type());
  PyModule_AddObject(module, "HFTransform", (PyObject *)LaDa::crystal::hftransform_type());
}
