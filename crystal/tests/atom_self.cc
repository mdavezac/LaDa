#include "LaDaConfig.h"

#include <python/numpy_types.h>
#include <boost/exception/get_error_info.hpp>
#include <boost/exception/diagnostic_information.hpp>

#include "../atom/atom.h"

namespace bp = boost::python;
using namespace LaDa::crystal;
static Atom satom; 
PyObject* get_static_object() { return satom.new_ref(); }
PyObject* set_static_object(PyObject* _module, PyObject *_object)
{
  try { satom.reset(_object); }
  catch(LaDa::error::pyerror &_e) { return NULL; }

  Py_RETURN_NONE; 
}

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

#ifdef LADA_DECLARE
#  error LADA_DECLARE already defined.
#endif
#define LADA_DECLARE(name, args) {#name, (PyCFunction)name, METH_ ## args, ""} 

static PyMethodDef methods[] = { 
  LADA_DECLARE(get_static_object, NOARGS),
  LADA_DECLARE(set_static_object, O),
  {NULL},
};

#undef LADA_DECLARE

PyMODINIT_FUNC init_atom_self(void) 
{
  import_array();
  LaDa::crystal::import();
  PyObject* module = Py_InitModule("_atom_self", methods);
}
