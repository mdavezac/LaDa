#include "LaDaConfig.h"

#include "../crystal.h"

using namespace LaDa::crystal;
static PyAtomObject *pysatom;
PyObject* get_static_object(PyObject* _module, PyObject*)
{ 
  std::cout << "static " << pysatom << std::endl;
  Py_INCREF(pysatom);
  return (PyObject*)pysatom;
}
PyObject* set_static_object(PyObject* _module, PyObject *_object)
{
  std::cout << "AM HERE 0 " << _module << std::endl;
  pysatom = NULL;
  std::cout << "AM HERE 1"  << std::endl;
  Atom satom = (PyAtomObject*)PyObject_GetAttrString(_module, "_atom");
  std::cout << "AM HERE 2"  << std::endl;
  try { satom.reset(_object); }
  catch(...) { std::cout << "Caught error" << std::endl; return NULL; }
  std::cout << "AM HERE 3"  << std::endl;
  pysatom = (PyAtomObject*)satom.borrowed();
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
  if(not LaDa::crystal::import()) return;
  Atom satom; 
  pysatom = (PyAtomObject*)satom.borrowed();
  std::cout << "Not null " << satom << std::endl;
  PyModule_AddObject(module, "_atom", (PyObject *)satom.new_ref());
}
