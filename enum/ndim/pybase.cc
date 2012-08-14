#include "LaDaConfig.h"

#include <Python.h>
#include <structmember.h>
#define PY_ARRAY_UNIQUE_SYMBOL enumeration_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>


#include <python/exceptions.h>
#include <python/wrap_numpy.h>
#include <python/quantity.h>

#include "pybase.h"
// NDimIterator member functions.
#include "cdi.hpp"

namespace LaDa
{
  namespace enumeration
  {
    //! Creates a new ndimiterator.
    NDimIterator* PyNDimIterator_New()
    {
      NDimIterator* result = (NDimIterator*) ndimiterator_type()->tp_alloc(ndimiterator_type(), 0);
      if(not result) return NULL;
      result->yielded = NULL;
      new(&result->ends) std::vector<t_ndim>;
      new(&result->vector) std::vector<t_ndim>;
      return result;
    }
    //! Creates a new ndimiterator with a given type.
    NDimIterator* PyNDimIterator_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs)
    {
      NDimIterator* result = (NDimIterator*)_type->tp_alloc(_type, 0);
      if(not result) return NULL;
      result->yielded = NULL;
      new(&result->ranges) std::vector< std::pair<t_ndim, t_ndim> >;
      new(&result->vector) std::vector<t_ndim>;
      return result;
    }

    // Creates a new ndimiterator with a given type, also calling initialization.
    NDimIterator* PyNDimIterator_NewFromArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs)
    {
      NDimIterator* result = PyNDimIterator_NewWithArgs(_type, _args, _kwargs);
      if(result == NULL) return NULL;
      if(_type->tp_init((PyObject*)result, _args, _kwargs) < 0) {Py_DECREF(result); return NULL; }
      return result;
    }

    // Returns pointer to ndimiterator type.
    PyTypeObject* ndimiterator_type()
    {
      static PySequenceMethods ndimiterator_as_sequence = {
          (lenfunc)ndimiterator_size,                    /* sq_length */
          (binaryfunc)NULL,                           /* sq_concat */
          (ssizeargfunc)NULL,                         /* sq_repeat */
          (ssizeargfunc)ndimiterator_getitem,            /* sq_item */
          (ssizessizeargfunc)NULL,                    /* sq_slice */
          (ssizeobjargproc)ndimiterator_setitem,         /* sq_ass_item */
          (ssizessizeobjargproc)NULL,                 /* sq_ass_slice */
          (objobjproc)NULL,                           /* sq_contains */
          (binaryfunc)NULL,                           /* sq_inplace_concat */
          (ssizeargfunc)NULL,                         /* sq_inplace_repeat */
      };
  
      static PyTypeObject dummy = {
          PyObject_HEAD_INIT(NULL)
          0,                                 /*ob_size*/
          "lada.enum.cppwrappers.NDimIterator",   /*tp_name*/
          sizeof(NDimIterator),              /*tp_basicsize*/
          0,                                 /*tp_itemsize*/
          (destructor)ndimiterator_dealloc,  /*tp_dealloc*/
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
          Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC | Py_TPFLAGS_HAVE_ITER, /*tp_flags*/
          "Defines a ndimiterator.\n\n"         /*tp_doc*/
          (traverseproc)ndimiterator_traverse,  /* tp_traverse */
          (inquiry)ndimiterator_gcclear,        /* tp_clear */
          0,		                     /* tp_richcompare */
          offsetof(NDimIterator, weakreflist),   /* tp_weaklistoffset */
          (getiterfunc)ndimiteratoriterator_self,  /* tp_iter */
          (iternextfunc)ndimiteratoriterator_next, /* tp_iternext */
          0,		                     /* tp_iternext */
          0,                                 /* tp_methods */
          0,                                 /* tp_members */
          0,                                 /* tp_getset */
          0,                                 /* tp_base */
          0,                                 /* tp_dict */
          0,                                 /* tp_descr_get */
          0,                                 /* tp_descr_set */
          0,                                 /* tp_dictoffset */
          (initproc)ndimiterator_init,       /* tp_init */
          0,                                 /* tp_alloc */
          (newfunc)PyNDimIterator_NewWithArgs,  /* tp_new */
      };
      if(dummy.tp_getattro == 0) dummy.tp_getattro = PyObject_GenericGetAttr;
      if(dummy.tp_setattro == 0) dummy.tp_setattro = PyObject_GenericSetAttr;
      return &dummy;
    }

  } // namespace Crystal
} // namespace LaDa
