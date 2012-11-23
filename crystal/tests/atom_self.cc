#include "LaDaConfig.h"

#include "../crystal.h"

using namespace LaDa::crystal;
PyObject* get_static_object(PyObject* _module, PyObject*)
{ 
  LaDa::python::Object module = PyImport_ImportModule("_atom_self");
  if(not module) return NULL;
  PyObject *result = PyObject_GetAttrString(module.borrowed(), "_atom");
  return result;
}
PyObject* set_static_object(PyObject* _module, PyObject *_object)
{
  if(not check_atom(_object))
  {
    LADA_PYERROR(TypeError, "Wrong type.");
    return NULL;
  }
  LaDa::python::Object module = PyImport_ImportModule("_atom_self");
  if(not module) return NULL;
  if( PyObject_SetAttrString(module.borrowed(), "_atom", _object) < 0) return NULL;
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
  if(not module) return; 
  if(not LaDa::python::import()) return;
  if(not LaDa::crystal::import()) return;
  Atom satom; 
  PyModule_AddObject(module, "_atom", (PyObject *)satom.new_ref());
  satom.release();
}
