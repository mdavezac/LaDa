#include "LaDaConfig.h"

#include <Python.h>
#include <structmember.h>
#define PY_ARRAY_UNIQUE_SYMBOL lada_math_ARRAY_API
#include <numpy/arrayobject.h>

#include <cmath>
#include <iterator> 
#include <limits>
#include <list>
#include <algorithm>

#include <boost/bind.hpp>
#include <boost/ref.hpp>

#include <Eigen/LU> 

#include <math/misc.h>
#include <math/fuzzy.h>
#include <math/eigen.h>
#include <math/smith_normal_form.h>
#include <math/quantity.h>
#include <math/gruber.h>

#include <errors/exceptions.h>
#include <python/python.h>

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif
#define LADA_CRYSTAL_MODULE 0
#include "crystal.h"

namespace LaDa
{
  namespace crystal
  {
    namespace 
    {
#     include "atom/pybase.cc"
#     include "structure/pybase.cc"
#     include "hart-forcade/pybase.cc"
#     include "utilities.cc"
#     include "map_sites.cc"
#     include "equivalent_structures.cc"
#     include "primitive.cc"
#     include "space_group.cc"
#     include "neighbors.cc"
#     include "coordination_shells.cc"
#     include "confsplit.cc"
#     include "periodic_dnc.cc"

#     include "methods.hpp"
    }
  }
}



PyMODINIT_FUNC initcppwrappers(void) 
{
  using namespace LaDa::crystal;
  static void *api_capsule[LADA_SLOT(crystal)];
  PyObject *c_api_object;

  char const doc[] =  "Wrapper around C++ atom/structure class and affiliates.";
  PyObject* module = Py_InitModule3("cppwrappers", methods_table, doc);
  if(not module) return;
  if(not LaDa::python::import()) return;
  import_array();

  /* Initialize the C API pointer array */
# undef PYLADA_CRYSTALMODULE_H
# include "crystal.h"

  /* Create a Capsule containing the API pointer array's address */
  static const char name[] = "lada.crystal.cppwrappers._C_API";
  c_api_object = PyCapsule_New((void *)api_capsule, name, NULL);
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
