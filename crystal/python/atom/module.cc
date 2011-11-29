#include "LaDaConfig.h"

#include <iterator> 

#include <Python.h>
#include <boost/python/tuple.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/borrowed.hpp>
#include <boost/python/slice.hpp>
#include <boost/python/slice_nil.hpp>

#include <math/python/python.hpp>
#include <python/exceptions.h>
#include <python/numpy_types.h>

#include "../../atom.h"
#include "extract_species.h"

namespace bp = boost::python;
namespace lp = LaDa::python;

using namespace LaDa;
using namespace LaDa::python;

// sequence interface
#include "sequence.hpp"

//! Structure holding shared pointer to an atom.
extern "C" struct AtomStr
{
  PyObject_HEAD
  PyArrayObject *position;
  PyObject* dictionary;
  PyObject* weakreflist;
  boost::shared_ptr< LaDa::crystal::AtomData< std::string > > atom;
};
//! Structure holding shared pointer to an atom.
extern "C" struct AtomSequence
{
  PyObject_HEAD
  PyArrayObject *position;
  Sequence *sequence;
  PyObject* dictionary;
  PyObject* weakreflist;
  boost::shared_ptr< LaDa::crystal::AtomData< std::vector<std::string> > > atom;
};


#ifdef LADA_NAME
#  error LADA_NAME already defined
#endif
#ifdef LADA_TYPE
#  error LADA_TYPE already defined
#endif
#ifdef LADA_ATOM_NUMBER
#  error LADA_ATOM_NUMBER already defined
#endif
#define LADA_ATOM_NUMBER 0
#define LADA_TYPE AtomStr
#define LADA_NAME(name) atomstr_ ## name

// atom getters/settters.
#include "getset.hpp"
// creation, deallocation, initialization.
#include "cdi.hpp"
// atom type declaration.
#include "type.hpp"

#undef LADA_NAME
#undef LADA_TYPE
#undef LADA_ATOM_NUMBER

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

static PyMethodDef atomstr_methods[] = 
  {
    {"_new_sequence", (PyCFunction)sequence_new_sequence, METH_VARARGS,
     "Creates new sequence. For debugging only.\n\n"
     "Sequences should only ever be obtained from an atom structure, and never created from scratch." },
    NULL
  };

PyMODINIT_FUNC initatom(void) 
{
  import_array(); // needed for NumPy 

  // set generic attributes for atom.
  atomstr_type.tp_getattro = PyObject_GenericGetAttr;
  atomstr_type.tp_setattro = PyObject_GenericSetAttr;

  if (PyType_Ready(&atomstr_type) < 0) return;
  if (PyType_Ready(&sequence_type) < 0) return;
  if (PyType_Ready(&sequenceiterator_type) < 0) return;

  Py_INCREF(&atomstr_type);
  Py_INCREF(&sequence_type);
  Py_INCREF(&sequenceiterator_type);

  char const doc[] =  "Wrapper around C++ atom class, and affiliates.";
  PyObject* module = Py_InitModule3("atom", atomstr_methods, doc);

  PyModule_AddObject(module, "AtomStr", (PyObject *)&atomstr_type);
  PyModule_AddObject(module, "Sequence", (PyObject *)&sequence_type);
  PyModule_AddObject(module, "SequenceIterator", (PyObject *)&sequenceiterator_type);
}
