#include "LaDaConfig.h"

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL lada_vff_ARRAY_API
#include <numpy/arrayobject.h>

#include <algorithm>

#include <errors/exceptions.h>
#include <crystal/crystal.h>
#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

#include "node/pybase.h"
#include "edge/pybase.h"
namespace LaDa
{
  namespace vff
  {
#   include "vff/zb.cc"
  }
}

namespace LaDa
{
  namespace vff
  {
    //! Methods table for vff module.
    static PyMethodDef methods_table[] = {
        { "_zb_energy", (PyCFunction) zincblende::energy, METH_VARARGS,
          "Computes energy for zinc-blende\n\n"
          "This is a very specialized function to compute energies for\n"
          "zinc-blende. It expects 6 bond-parameters and 7 angle-parameters in\n"
          "numpy arrays of type 'float64'. Anything else may have unpredictable\n"
          "results." },
        { "_zb_jacobian", (PyCFunction) zincblende::jacobian, METH_VARARGS,
          "Computes jacobian for zinc-blende\n\n"
          "This is a very specialized function to compute jacobian for\n"
          "zinc-blende. It expects 6 bond-parameters and 7 angle-parameters in\n"
          "numpy arrays of type 'float64'. Anything else may have unpredictable\n"
          "results." },
        {NULL, NULL, 0, NULL}        /* Sentinel */
    };
  }
}
PyMODINIT_FUNC initcppwrappers(void) 
{
  char const doc[] =  "Wrapper around C++ vff class and affiliates.";
  PyObject* module = Py_InitModule3("cppwrappers", LaDa::vff::methods_table, doc);
  if(not module) return;
  import_array(); // needed for NumPy 
  if(not LaDa::python::import()) return;
  if(not LaDa::math::import()) return;
  if(not LaDa::crystal::import()) return;

  if (PyType_Ready(LaDa::vff::node_type()) < 0) return;
  if (PyType_Ready(LaDa::vff::edge_type()) < 0) return;
  if (PyType_Ready(LaDa::vff::bonditerator_type()) < 0) return;
  if (PyType_Ready(LaDa::vff::dcbonditerator_type()) < 0) return;
  if (PyType_Ready(LaDa::vff::angleiterator_type()) < 0) return;

  Py_INCREF(LaDa::vff::node_type());
  Py_INCREF(LaDa::vff::edge_type());
  Py_INCREF(LaDa::vff::bonditerator_type());
  Py_INCREF(LaDa::vff::dcbonditerator_type());
  Py_INCREF(LaDa::vff::angleiterator_type());


  PyModule_AddObject(module, "Node", (PyObject *)LaDa::vff::node_type());
  PyModule_AddObject(module, "Edge", (PyObject *)LaDa::vff::edge_type());
  PyModule_AddObject(module, "BondIterator", (PyObject *)LaDa::vff::bonditerator_type());
  PyModule_AddObject(module, "ScBondIterator", (PyObject *)LaDa::vff::dcbonditerator_type());
  PyModule_AddObject(module, "AngleIterator", (PyObject *)LaDa::vff::angleiterator_type());
}
