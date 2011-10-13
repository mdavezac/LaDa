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

#include "../atom.h"
#include "extract_species.h"

namespace bp = boost::python;
namespace lp = LaDa::python;

using namespace LaDa;
using namespace LaDa::python;

//! Structure holding shared pointer to an atom.
extern "C" struct AtomStr
{
  PyObject_HEAD
  boost::shared_ptr< LaDa::crystal::AtomData< std::string > > atom;
};


// atom getters/settters.
#include "atom_getset.hpp"
// creation, deallocation, initialization.
#include "atom_cdi.hpp"
// atom type declaration.
#include "atom_type.hpp"

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

static PyMethodDef AtomStr_methods[] = { {NULL} };

PyMODINIT_FUNC initatom(void) 
{
  import_array(); // needed for NumPy 
  PyObject* m;

  if (PyType_Ready(&AtomStrType) < 0) return;

  m = Py_InitModule3("atom", AtomStr_methods,
                     "Example module that creates an extension type.");

  Py_INCREF(&AtomStrType);
  PyModule_AddObject(m, "Atom", (PyObject *)&AtomStrType);
}
