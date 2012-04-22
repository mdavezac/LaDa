#include "LaDaConfig.h"

#include <Python.h>
#include <numpy/arrayobject.h>
#include <boost/exception/get_error_info.hpp>
#include <boost/exception/diagnostic_information.hpp>

#include <opt/debug.h>

#include "../quantity.h"


namespace bp = boost::python;
using namespace LaDa::python;
PyObject* get_static_object()
{ 
  try { return Quantity(1, "m").new_ref(); }
  catch(...) { return NULL; }
}
PyObject* check_get() 
{
  LADA_DOASSERT( std::abs(Quantity(1, "m").get("cm") -1e-2) < 1e-12, "Did not get correct scale."); 
  LADA_DOASSERT( std::abs(Quantity(1, "m").get() -1e0) < 1e-12, "Did not get correct scale."); 
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
  LADA_DECLARE(check_get, NOARGS),
  {NULL},
};

#undef LADA_DECLARE

PyMODINIT_FUNC init_quantity(void) 
{
  import_array();
  {
    Object o(PyImport_ImportModule("lada.error")); 
    if(not o) return;
  }
  PyObject* module = Py_InitModule("_quantity", methods);
}