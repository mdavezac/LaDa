#include "LaDaConfig.h"

#include <Python.h>
#include <structmember.h>
#define PY_ARRAY_UNIQUE_SYMBOL lada_math_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>


#include <python/exceptions.h>
#include <crystal/python/wrap_numpy.h>
#include <math/fuzzy.h>
#include <math/smith_normal_form.h>

#include "../structure/pybase.h"

#include "pybase.h"
// hftransform getters/settters.
#include "get.hpp"
// hftransform member functions.
#include "members.hpp"
// creation, deallocation, initialization.
#include "cdi.hpp"

namespace LaDa
{
  namespace crystal
  {
    namespace // limits declarations to current namespace.
    {
      //! Creates a new hftransform with a given type.
      PyHFTObject* PyHFTransform_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs)
      {
        PyHFTObject* result = (PyHFTObject*)_type->tp_alloc(_type, 0);
        return result;
      }
      
      // Creates a new hftransform with a given type, also calling initialization.
      PyHFTObject* new_hftransform(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs)
      {
        PyHFTObject* result = PyHFTransform_NewWithArgs(_type, _args, _kwargs);
        if(result == NULL) return NULL;
        if(_type->tp_init((PyObject*)result, _args, _kwargs) < 0) {Py_DECREF(result); return NULL; }
        return result;
      }
      
      // Creates a deepcopy of hftransform.
      PyHFTObject *copy_hftransform(PyHFTObject* _self, PyObject *_memo)
      {
        PyHFTObject* result = (PyHFTObject*)_self->ob_type->tp_alloc(_self->ob_type, 0);
        if(not result) return NULL;
        new(&result->transform) LaDa::math::rMatrix3d(_self->transform);
        new(&result->quotient) LaDa::math::iVector3d(_self->quotient);
        return result;
      }
      // Returns pointer to hftransform type.
      PyTypeObject* hftransform_type()
      {
#       ifdef LADA_DECLARE
#         error LADA_DECLARE already defined.
#       endif
#       define LADA_DECLARE(name, doc) \
          {const_cast<char*>(#name), (getter) hftransform_get ## name, NULL, const_cast<char*>(doc)}
        
        static PyGetSetDef getsetters[] = {
            LADA_DECLARE(transform, "Transformation matrix to go from position to hf index."),
            LADA_DECLARE(quotient, "Periodicity quotient."),
            LADA_DECLARE(size, "Number of unitcells in the supercell."),
            {NULL}  /* Sentinel */
        };
#       undef LADA_DECLARE
#       define LADA_DECLARE(name, func, args, doc) \
          {#name, (PyCFunction)func, METH_ ## args, doc} 
        static PyMethodDef methods[] = {
            LADA_DECLARE(copy, hftransform_copy, NOARGS, "Returns a deepcopy of the hftransform."),
            LADA_DECLARE(__copy__, hftransform_shallowcopy, NOARGS, "Shallow copy of an hftransform."),
            LADA_DECLARE(__deepcopy__, copy_hftransform, O, "Deep copy of an hftransform."),
            LADA_DECLARE(__getstate__, hftransform_getstate, NOARGS, "Implements pickle protocol."),
            LADA_DECLARE(__setstate__, hftransform_setstate, O, "Implements pickle protocol."),
            LADA_DECLARE(__reduce__,   hftransform_reduce, NOARGS, "Implements pickle protocol."),
            LADA_DECLARE( indices,   hftransform_indices, O,
                          "indices of input atomic position in cyclic Z-group.\n\n"
                          ":param vec: A 3d-vector in the sublattice of interest.\n"
                          ":returns: The 3-d indices in the cyclic group.\n" ),
            LADA_DECLARE( flatten_indices,   hftransform_flatten_indices, VARARGS | METH_KEYWORDS,
                          "Flattens cyclic Z-group indices.\n\n"
                          ":param int i: \n"
                          "  First index into cyclic Z-group.\n"
                          ":param int j: \n"
                          "  Second index into cyclic Z-group.\n"
                          ":param int j: \n"
                          "  Third index into cyclic Z-group.\n"
                          ":param int site: \n"
                          "  Optional site index for multilattices.\n"
                          ":returns: An integer which can serve as an index into a 1d array.\n"), 
            LADA_DECLARE( index,   hftransform_flat_index, VARARGS | METH_KEYWORDS,
                          "Flat index into cyclic Z-group.\n\n"
                          ":param pos: (3d-vector)\n"
                          "    Atomic position with respect to the sublattice of\n"
                          "    interest.  Do not forget to shift the sublattice\n"
                          "    back to the origin.\n"
                          ":param site: (integer)\n"
                          "   Optional site index. If there are more than one\n"
                          "   sublattice in the structure, then the flattened\n"
                          "   indices need to take this into account.\n" ),
            {NULL}  /* Sentinel */
        };
#       undef LADA_DECLARE
      
        static PyTypeObject dummy = {
            PyObject_HEAD_INIT(NULL)
            0,                                 /*ob_size*/
            "lada.crystal.cppwrappers.HFTransform",   /*tp_name*/
            sizeof(PyHFTObject),             /*tp_basicsize*/
            0,                                 /*tp_itemsize*/
            0,                                 /*tp_dealloc*/
            0,                                 /*tp_print*/
            0,                                 /*tp_getattr*/
            0,                                 /*tp_setattr*/
            0,                                 /*tp_compare*/
            0,                                 /*tp_repr*/
            0,                                 /*tp_as_number*/
            0,                                 /*tp_as_sequence*/
            0,                                 /*tp_as_mapping*/
            0,                                 /*tp_hash */
            0,                                 /*tp_call*/
            0,                                 /*tp_str*/
            0,                                 /*tp_getattro*/
            0,                                 /*tp_setattro*/
            0,                                 /*tp_as_buffer*/
            Py_TPFLAGS_DEFAULT,                /*tp_flags*/
            "Defines a hftransform.\n\n"    /*tp_doc*/
            "The Hart-Forcade transform computes the cyclic group of supercell\n"
            "with respect to its backbone lattice. It can then be used to index\n"
            "atoms in the supercell, irrespective of which periodic image is given\n"
            "[HF]_.\n\n"
            ":param lattice:\n"
            "   Defines the cyclic group.\n"
            ":type lattice: :py:func:`Structure` or matrix\n"
            ":param supercell:\n"
            "    Supercell for which to compute cyclic group.\n"
            ":type supercell: :py:func:`Structure` or matrix\n",
            0,                                 /* tp_traverse */
            0,                                 /* tp_clear */
            0,		                     /* tp_richcompare */
            0,                                 /* tp_weaklistoffset */
            0,                                 /* tp_iter */
            0,		                     /* tp_iternext */
            methods,                           /* tp_methods */
            0,                                 /* tp_members */
            getsetters,                        /* tp_getset */
            0,                                 /* tp_base */
            0,                                 /* tp_dict */
            0,                                 /* tp_descr_get */
            0,                                 /* tp_descr_set */
            0,                                 /* tp_dictoffset */
            (initproc)hftransform_init,     /* tp_init */
            0,                                 /* tp_alloc */
            (newfunc)PyHFTransform_NewWithArgs,  /* tp_new */
        };
        return &dummy;
      }
    } // anonymous namespace
  } // namespace crystal
} // namespace LaDa
