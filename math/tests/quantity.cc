#include "LaDaConfig.h"

#include "../quantity.h"


using namespace LaDa::math;
PyObject* is_quantity(PyObject *_module, PyObject *_in)
{ 
  if(not PyQuantity_Check(_in)) Py_RETURN_FALSE;
  Py_RETURN_TRUE;
}
PyObject* fromC(PyObject *_module, PyObject *_args)
{
  double pynumber;
  char *pyunits;
  if(not PyArg_ParseTuple(_args, "ds", &pynumber, &pyunits)) return NULL;
  return PyQuantity_FromC(pynumber, std::string(pyunits));
}
PyObject* fromPy(PyObject *_module, PyObject *_args)
{
  PyObject *number;
  PyObject *units;
  if(not PyArg_ParseTuple(_args, "OO", &number, &units)) return NULL;
  return PyQuantity_FromPy(number, units);
}
PyObject* get_angstrom(PyObject *_module, PyObject *_in)
{
  LaDa::types::t_real result(PyQuantity_Get(_in, "angstrom"));
  if(std::abs(result) < 1e-8 and PyErr_Occurred())
  {
    PyErr_Clear();
    Py_RETURN_NONE;
  }
  return PyFloat_FromDouble(result);
}

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

#ifdef LADA_DECLARE
#  error LADA_DECLARE already defined.
#endif
#define LADA_DECLARE(name, args) {#name, (PyCFunction)name, METH_ ## args, ""} 

static PyMethodDef methods[] = { 
  LADA_DECLARE(is_quantity, O),
  LADA_DECLARE(fromC, VARARGS),
  LADA_DECLARE(fromPy, VARARGS),
  LADA_DECLARE(get_angstrom, O),
  {NULL},
};

#undef LADA_DECLARE

PyMODINIT_FUNC init_quantity(void) 
{
  PyObject* module = Py_InitModule("_quantity", methods);
}
