#include "LaDaConfig.h"

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

//! Structure holding shared pointer to an atom.
extern "C" struct AtomStr
{
  PyObject_HEAD
  PyArrayObject *position;
  PyObject* dictionary;
  PyObject* weakreflist;
  boost::shared_ptr< LaDa::crystal::AtomData< std::string > > atom;
};


// atom getters/settters.
#include "getset.hpp"
// creation, deallocation, initialization.
#include "cdi.hpp"
// atom type declaration.
#include "type.hpp"
// set interface
#include "set.hpp"

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

static PyMethodDef atomstr_methods[] = 
  {
    {"_new_set", (PyCFunction)set_new_set, METH_VARARGS,
     "Creates new set. For debugging only.\n\n"
     "Sets should only ever be obtained from an atom structure, and never created from scratch." },
    NULL
  };

PyMODINIT_FUNC initatom(void) 
{
  import_array(); // needed for NumPy 

  atomstr_type.tp_getattro = PyObject_GenericGetAttr;
  atomstr_type.tp_setattro = PyObject_GenericSetAttr;

  if (PyType_Ready(&atomstr_type) < 0) return;
  if (PyType_Ready(&set_type) < 0) return;
  if (PyType_Ready(&setiterator_type) < 0) return;

  Py_INCREF(&atomstr_type);
  Py_INCREF(&set_type);
  Py_INCREF(&setiterator_type);

  char const doc[] =  "Wrapper around C++ atom class, and affiliates.";
  PyObject* module = Py_InitModule3("atom", atomstr_methods, doc);

  PyModule_AddObject(module, "AtomStr", (PyObject *)&atomstr_type);
  PyModule_AddObject(module, "Set", (PyObject *)&set_type);
  PyModule_AddObject(module, "SetIterator", (PyObject *)&setiterator_type);
}
