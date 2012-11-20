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
#define LADA_CRYSTAL_MODULE
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
  static void *api_capsule[26];
  PyObject *c_api_object;

  char const doc[] =  "Wrapper around C++ atom/structure class and affiliates.";
  PyObject* module = Py_InitModule3("cppwrappers", methods_table, doc);
  if(not module) return;

  import_array(); // needed for NumPy 

  std::cout << BOOST_PP_SLOT(1) << std::endl;
  /* Initialize the C API pointer array */
  api_capsule[0] = (void *)atom_type();
  api_capsule[1] = (void *)((PyAtomObject*(*)())new_atom);
  api_capsule[2] = (void *)((PyAtomObject*(*)(PyTypeObject*, PyObject*, PyObject*))new_atom);
  api_capsule[3] = (void *)copy_atom;

  api_capsule[4] = (void *)structure_type();
  api_capsule[5] = (void *)structureiterator_type();
  api_capsule[6] = (void *)((PyStructureObject*(*)())new_structure);
  api_capsule[7] = (void *)((PyStructureObject*(*)(PyTypeObject*, PyObject*, PyObject*))new_structure);
  api_capsule[8] = (void *)copy_structure;
  api_capsule[9] = (void *)itransform_structure;

  api_capsule[20] = (void *)hftransform_type();
  api_capsule[21] = (void *)new_hftransform;
  api_capsule[22] = (void *)copy_hftransform;
  api_capsule[23] = (void *)_init_hft;

  api_capsule[10] = (void *)map_sites;
  api_capsule[11] = (void *)equivalent;
  api_capsule[12] = (void *)primitive;
  api_capsule[13] = (void *)is_primitive;
  api_capsule[14] = (void *)cell_invariants;
  api_capsule[15] = (void *)space_group;
  api_capsule[16] = (void *)neighbors;
  api_capsule[17] = (void *)coordination_shells;
  api_capsule[18] = (void *)splitconfigs;
  api_capsule[19] = (void *)dnc_boxes;

# ifdef LADA_TYPEDEF
#  error LADA_TYPEDEF already defined
# endif
# define LADA_TYPEDEF                                                          \
     (LaDa::math::rVector3d(*)( LaDa::math::rVector3d const&,                  \
                                LaDa::math::rMatrix3d const&,                  \
                                LaDa::math::rMatrix3d const& ))
  api_capsule[24] = (void *)(LADA_TYPEDEF into_cell);
  api_capsule[25] = (void *)(LADA_TYPEDEF into_voronoi);
  api_capsule[26] = (void *)(LADA_TYPEDEF zero_centered);
# undef LADA_TYPEDEF
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
