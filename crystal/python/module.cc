#include "LaDaConfig.h"

#include <iterator> 

#include <Python.h>

#include <math/python/python.hpp>
#include <python/exceptions.h>
#include <python/numpy_types.h>

#include "../cast.h"

#include "atom/atom.h"
#include "structure/structure.h"

namespace lp = LaDa::python;

using namespace LaDa;
using namespace LaDa::python;

// sequence interface
#include "atom/sequence.hpp"

#ifdef LADA_ATOM_CTYPE
#  error LADA_ATOM_CTYPE already defined
#endif
#ifdef LADA_CTYPE
#  error LADA_CTYPE already defined
#endif
#ifdef LADA_NAME
#  error LADA_NAME already defined
#endif
#ifdef LADA_TYPE
#  error LADA_TYPE already defined
#endif
#ifdef LADA_ATOM_NAME
#  error LADA_NAME already defined
#endif
#ifdef LADA_ATOM_TYPE
#  error LADA_TYPE already defined
#endif
#ifdef LADA_ATOM_NUMBER
#  error LADA_ATOM_NUMBER already defined
#endif
#define LADA_ATOM_NUMBER 0
#define LADA_TYPE StructureStr
#define LADA_CTYPE LaDa::crystal::StructureData<std::string>
#define LADA_NAME(name) structurestr_ ## name
#define LADA_ATOM_TYPE AtomStr
#define LADA_ATOM_CTYPE LaDa::crystal::AtomData<std::string> 
#define LADA_ATOM_NAME(name) atomstr_ ## name
// atom getters/settters.
#include "atom/getset.hpp"
// atom member functions.
#include "atom/members.h"
// creation, deallocation, initialization.
#include "atom/cdi.hpp"
// atom type declaration.
#include "atom/type.hpp"
// structure getters/settters.
#include "structure/getset.h"
// structure member functions.
#include "structure/members.h"
// creation, deallocation, initialization.
#include "structure/cdi.h"
// Mapping methods.
#include "structure/mapping.h"
// structure type declaration.
#include "structure/type.h"
#undef LADA_ATOM_NAME
#undef LADA_ATOM_TYPE
#undef LADA_ATOM_CTYPE
#undef LADA_NAME
#undef LADA_TYPE
#undef LADA_CTYPE
#undef LADA_ATOM_NUMBER
#define LADA_ATOM_NUMBER 1
#define LADA_ATOM_CTYPE LaDa::crystal::AtomData< std::vector<std::string> >
#define LADA_ATOM_TYPE AtomSequence
#define LADA_ATOM_NAME(name) atomsequence_ ## name
#define LADA_CTYPE LaDa::crystal::StructureData< std::vector<std::string> >
#define LADA_TYPE StructureSequence
#define LADA_NAME(name) structuresequence_ ## name
// atom getters/settters.
#include "atom/getset.hpp"
// atom member functions.
#include "atom/members.h"
// creation, deallocation, initialization.
#include "atom/cdi.hpp"
// atom type declaration.
#include "atom/type.hpp"
// structure getters/settters.
#include "structure/getset.h"
// structure member functions.
#include "structure/members.h"
// creation, deallocation, initialization.
#include "structure/cdi.h"
// Mapping methods.
#include "structure/mapping.h"
// structure type declaration.
#include "structure/type.h"
#undef LADA_NAME
#undef LADA_ATOM_TYPE
#undef LADA_ATOM_CTYPE
#undef LADA_ATOM_NUMBER
#undef LADA_TYPE
#undef LADA_CTYPE
#undef LADA_NUMBER

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

static PyMethodDef crystal_methods[] = 
  {
    {"_new_sequence", (PyCFunction)sequence_new_sequence, METH_VARARGS,
     "Creates new sequence. For debugging only.\n\n"
     "Sequences should only ever be obtained from an atom structure, and never created from scratch." },
    NULL
  };

PyMODINIT_FUNC initcppwrappers(void) 
{
  import_array(); // needed for NumPy 

  if (PyType_Ready(&structurestr_type) < 0) return;
  if (PyType_Ready(&structuresequence_type) < 0) return;
  if (PyType_Ready(&atomstr_type) < 0) return;
  if (PyType_Ready(&atomsequence_type) < 0) return;
  if (PyType_Ready(&sequence_type) < 0) return;
  if (PyType_Ready(&sequenceiterator_type) < 0) return;

  Py_INCREF(&structurestr_type);
  Py_INCREF(&structuresequence_type);
  Py_INCREF(&atomstr_type);
  Py_INCREF(&atomsequence_type);
  Py_INCREF(&sequence_type);
  Py_INCREF(&sequenceiterator_type);

  char const doc[] =  "Wrapper around C++ atom/structure class and affiliates.";
  PyObject* module = Py_InitModule3("cppwrappers", crystal_methods, doc);

  PyModule_AddObject(module, "StructureStr", (PyObject *)&structurestr_type);
  PyModule_AddObject(module, "StructureSequence", (PyObject *)&structuresequence_type);
  PyModule_AddObject(module, "AtomStr", (PyObject *)&atomstr_type);
  PyModule_AddObject(module, "AtomSequence", (PyObject *)&atomsequence_type);
  PyModule_AddObject(module, "Sequence", (PyObject *)&sequence_type);
  PyModule_AddObject(module, "SequenceIterator", (PyObject *)&sequenceiterator_type);
}
