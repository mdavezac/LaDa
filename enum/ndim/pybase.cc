#include "LaDaConfig.h"

#include <Python.h>
#include <structmember.h>
#define PY_ARRAY_UNIQUE_SYMBOL enumeration_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>

#include <vector>
#include <algorithm>


#include <python/numpy_types.h>
#include <python/exceptions.h>
#include <python/wrap_numpy.h>
#include <python/quantity.h>

#include "pybase.h"

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

    // Function to deallocate a string atom.
    static void ndimiterator_dealloc(NDimIterator *_self)
    {
      gcclear(_self);
  
      // calls destructor explicitely.
      PyTypeObject* ob_type = _self->ob_type;
      _self->~NDimIterator();

      ob_type->tp_free((PyObject*)_self);
    }
  
    // Function to initialize an atom.
    static int ndimiterator_init(NDimIterator* _self, PyObject* _args, PyObject *_kwargs)
    {
      if(_args == NULL)
      {
        LADA_PYERROR(TypeError, "NDimIterator expects at least one argument.");
        return -1;
      }
      Py_ssize_t const N = PyTuple_Size(_args);
      if(N == 0)
      {
        LADA_PYERROR(TypeError, "NDimIterator expects at least one argument.");
        return -1;
      }
      if(_kwargs != NULL and PyDict_Size(_kwargs))
      {
        LADA_PYERROR(TypeError, "NDimIterator does not expect keyword arguments.");
        return -1;
      }
      _self.ends.clear();
      _self.ends.reserve(N);
      for(Py_ssize_t i(0); i < N; ++i)
      {
        PyObject *item = PyTuple_GET_ITEM(_args, i);
        if(PyLong_Check(item)) _self->ends.push_back( PyLong_AS_LONG(item));
        else if(PyInt_Check(item)) _self->ends.push_back( PyInt_AS_LONG(item));
        else
        {
          LADA_PYERROR(TypeError, "Unknown type in NDimIterator.");
          return -1;
        }
        if(_self.ends.back() < 0)
        {
          LADA_PYERROR(TypeError, "NDimIterator does not expect negative arguments.");
          return -1;
        }
        if(_self.ends.back() > 5)
        {
          LADA_PYERROR(TypeError, "Argument larger than 5. "
                       "Recompile LaDa if that is what you want.");
          return -1;
        }
      } 
      _self->vector.resize(_self->ends.size());
      std::fill(_self->vector.begin(), _self->vector.end(), 1);
      if(_self->yielded != NULL)
      {
        PyObject* dummy = _self.yielded;
        _self.yielded = NULL;
        Py_DECREF(dummy);
      }
      typedef math::numpy::type<t_ndim> t_type;
      _self->yielded = (PyArrayObject*) PyArray_SimpleNewFromData(1, &N, t_type::value, &counter[0]);
      if(not _self.yielded) return -1;
#     ifdef LADA_MACRO
#       error LADA_MACRO already defined
#     endif
#     ifdef NPY_ARRAY_WRITEABLE
#       define LADA_MACRO NPY_ARRAY_WRITEABLE
#     else
#       define LADA_MACRO NPY_WRITEABLE
#     endif
      if(_self->yielded->flags & LADA_MACRO) _self->yielded->flags -= LADA_MACRO;
#     undef LADA_MACRO
#     ifdef NPY_ARRAY_C_CONTIGUOUS
#       define LADA_MACRO NPY_ARRAY_C_CONTIGUOUS;
#     else 
#       define LADA_MACRO NPY_C_CONTIGUOUS
#     endif
#     undef LADA_MACRO
      if(not (_self->yielded->flags & LADA_MACRO)) _self->yielded->flags += LADA_MACRO;
      _self->yielded->parent = _self;
      Py_INCREF(_self);
      return 0;
    }
  
    static int traverse(NDimIterator *self, visitproc visit, void *arg)
    {
      Py_VISIT(self->yielded);
      return 0;
    }
  
    static int gcclear(NDimIterator *self)
    { 
      Py_CLEAR(self->yielded);
      return 0;
    }

    static PyObject* self(PyObject* _self)
    {
      Py_INCREF(_self);
      return _self;
    }

    static PyObject* iter(NDimIterator* _self)
    {
#     ifdef LADA_DEBUG
        if(_self->ends.size() != _self.counter.size())
        {
          LADA_PYTHROW(internal, "Counter and Ends have different size.");
          return NULL;
        }
#     endif
      std::vector<t_ndim>::iterator i_first = _self->counter.begin();
      std::vector<t_ndim>::iterator i_last = _self->counter.end();
      std::vector<t_ndim>::const_iterator i_end = _self->ends.begin();
      for(; i_first != i_last; ++i_first, ++i_end)
        if(*i_first == *i_end) *i_first = 1;
        else
        {
          ++(*i_first);
#         ifdef LADA_DEBUG
            if(_self->yielded == NULL)
            {
              LADA_PYTHROW(internal, "Yielded was not initialized.");
              return NULL;
            }
            if(_self->yielded->data != &_self.counter[0])
            {
              LADA_PYTHROW(internal, "Yielded does not reference counter.");
              return NULL;
            }
#         endif
          Py_INCREF(_self->yielded);
          return _self->yielded;
        }
      PyObject* dummy = (PyObject*) _self->yielded;
      _self->yielded = NULL;
      Py_DECREF(dummy);
      _self->ends.clear();
      _self->counter.clear();
      PyErr_SetNone(PyExc_StopIteration);
      return NULL;
    }


    // Returns pointer to ndimiterator type.
    PyTypeObject* ndimiterator_type()
    {
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
          "Defines a ndimiterator.\n\n"      /*tp_doc*/
          (traverseproc)traverse,            /* tp_traverse */
          (inquiry)gcclear,                  /* tp_clear */
          0,		                             /* tp_richcompare */
          0,                                 /* tp_weaklistoffset */
          (getiterfunc)self,                 /* tp_iter */
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
      return &dummy;
    }

  } // namespace Crystal
} // namespace LaDa
