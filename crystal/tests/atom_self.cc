#include "LaDaConfig.h"

#include "../atom/atom.h"

using namespace LaDa::crystal;
static Atom satom; 
PyObject* get_static_object() { return satom.new_ref(); }
PyObject* set_static_object(PyObject* _module, PyObject *_object)
{
  try { satom.reset(_object); }
  catch(...) { return NULL; }
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
  PyObject* module = Py_InitModule("_atom_self", methods);
}
