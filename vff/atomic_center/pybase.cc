#include "LaDaConfig.h"

#include <Python.h>
#include <structmember.h>
#define PY_ARRAY_UNIQUE_SYMBOL lada_vff_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>


#include <python/exceptions.h>
#include <python/wrap_numpy.h>
#include <math/quantity.h>

#include "pybase.h"
// atomic_center getters/settters.
#include "getset.hpp"
// iterator functions and type.
#include "iterator.hpp"
// atomic_center member functions.
#include "members.hpp"
// creation, deallocation, initialization.
#include "cdi.hpp"
// sequence functions.
#include "sequence.hpp"

namespace LaDa
{
  namespace vff
  {
    //! returns 1 angstrom.
    static PyObject* get_unit_angstrom();
    //! Creates a new atomic_center.
    AtomicCenterData* PyAtomicCenter_New()
    {
      AtomicCenterData* result = (AtomicCenterData*) atomic_center_type()->tp_alloc(atomic_center_type(), 0);
      if(not result) return NULL;
      result->weakreflist = NULL;
      new(&result->bonds) std::vector<AtomicCenterData*>;
      result->center = NULL;
      result->index = 0u;
      return result;
    }
    //! Creates a new atomic_center with a given type.
    AtomicCenterData* PyAtomicCenter_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs)
      { return PyAtomicCenter_New(); }

    // Creates a new atomic_center with a given type, also calling initialization.
    AtomicCenterData* PyAtomicCenter_NewFromArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs)
    {
      AtomicCenterData* result = PyAtomicCenter_NewWithArgs(_type, _args, _kwargs);
      if(result == NULL) return NULL;
      if(_type->tp_init((PyObject*)result, _args, _kwargs) < 0) {Py_DECREF(result); return NULL; }
      return result;
    }

    // Returns pointer to atomic_center type.
    PyTypeObject* atomic_center_type()
    {
      static PyMappingMethods atomic_center_as_mapping = {
          (lenfunc)atomic_center_size,
          (binaryfunc)atomic_center_subscript,
          (objobjargproc)atomic_center_ass_subscript
      };
      static PySequenceMethods atomic_center_as_sequence = {
          (lenfunc)atomic_center_size,                    /* sq_length */
          (binaryfunc)NULL,                           /* sq_concat */
          (ssizeargfunc)NULL,                         /* sq_repeat */
          (ssizeargfunc)atomic_center_getitem,            /* sq_item */
          (ssizessizeargfunc)NULL,                    /* sq_slice */
          (ssizeobjargproc)atomic_center_setitem,         /* sq_ass_item */
          (ssizessizeobjargproc)NULL,                 /* sq_ass_slice */
          (objobjproc)NULL,                           /* sq_contains */
          (binaryfunc)NULL,                           /* sq_inplace_concat */
          (ssizeargfunc)NULL,                         /* sq_inplace_repeat */
      };
#     ifdef LADA_DECLARE
#       error LADA_DECLARE already defined.
#     endif
#     define LADA_DECLARE(name, doc) \
        { const_cast<char*>(#name), (getter) atomic_center_get ## name, \
          (setter) atomic_center_set ## name, const_cast<char*>(doc) }
      
      static PyGetSetDef getsetters[] = {
          LADA_DECLARE(cell, "Cell matrix in cartesian coordinates.\n\n"
                             "Unlike most ab-initio codes, cell-vectors are "
                             "given in column vector format. " 
                             "The cell does not yet have units. "
                             "Units depend upon :class:`AtomicCenter.scale`. "
                             "Across lada, it is expected that a cell time "
                             "this scale are angstroms. Finally, the cell "
                             "is owned internally by the atomic_center. It cannot be set to "
                             "reference an object (say a list or numpy array). "
                             "``atomic_center.cell = some_list`` will copy the values of ``some_list``. "),
          LADA_DECLARE(scale, "Scale factor of this atomic_center.\n"
                              "Should be a number or unit given by "
                              "the python package `quantities "
                              "<http://packages.python.org/quantities/index.html>`_.\n\n"
                              ".. note:: The scale is always converted to angstroms. "),
          { const_cast<char*>("volume"), (getter) atomic_center_getvolume, NULL, 
            const_cast<char*>("Volume of the atomic_center.\n\nIncludes scale.") },
          {NULL}  /* Sentinel */
      };
#     undef LADA_DECLARE
#     define LADA_DECLARE(name, object, doc) \
        { const_cast<char*>(#name), T_OBJECT_EX, \
          offsetof(AtomicCenterData, object), 0, const_cast<char*>(doc) }
      static PyMemberDef members[] = {
        LADA_DECLARE(__dict__, pydict, "Python attribute dictionary."),
#       ifdef LADA_DEBUG
          LADA_DECLARE(_weakreflist, weakreflist, "List of weak references."),
#       endif
        {NULL, 0, 0, 0, NULL}  /* Sentinel */
      };
#     undef LADA_DECLARE
#     define LADA_DECLARE(name, func, args, doc) \
        {#name, (PyCFunction)func, METH_ ## args, doc} 
      static PyMethodDef methods[] = {
          LADA_DECLARE(copy, atomic_center_copy, NOARGS, "Returns a deepcopy of the atomic_center."),
          LADA_DECLARE( to_dict, atomic_center_to_dict, NOARGS, 
                        "Returns a dictionary with shallow copies of items." ),
          LADA_DECLARE(__copy__, atomic_center_shallowcopy, NOARGS, "Shallow copy of an atomic_center."),
          LADA_DECLARE(__deepcopy__, PyAtomicCenter_Copy, O, "Deep copy of an atomic_center."),
          LADA_DECLARE(__getstate__, atomic_center_getstate, NOARGS, "Implements pickle protocol."),
          LADA_DECLARE(__setstate__, atomic_center_setstate, O, "Implements pickle protocol."),
          LADA_DECLARE(__reduce__,   atomic_center_reduce, NOARGS, "Implements pickle protocol."),
          LADA_DECLARE( add_atom,   atomic_center_add_atom, KEYWORDS,
                        "Adds atom to atomic_center.\n\n"
                        "The argument to this function is either another atom, "
                        "in which case a reference to that atom is appended to "
                        "the atomic_center. Or, it is any arguments used to "
                        "initialize atoms in :class:`Atom`. "
                        "Finally, this function can be chained as follows::\n\n"
                        "  atomic_center.add_atom(0,0,0, 'Au')                        \\\n"
                        "           .add_atom(0.25, 0.25, 0.25, ['Pd', 'Si'], m=5)\\\n"
                        "           .add_atom(atom_from_another_atomic_center)\n\n"
                        "In the example above, both ``atomic_center`` and the *other* atomic_center will "
                        "reference the same atom (``atom_from_another_atomic_center``). "
                        "Changing, say, that atom's type in one atomic_center will also "
                        "change it in the other.\n\n"
                        ":returns: The atomic_center itself, so that add_atom methods can be chained."),
          LADA_DECLARE( insert, atomic_center_insert, VARARGS, 
                        "Inserts atom at given position.\n\n"
                        ":param index:\n    Position at which to insert the atom.\n"
                        ":param atom: \n    :class:`Atom` or subtype to insert." ),
          LADA_DECLARE(pop, atomic_center_pop, O, "Removes and returns atom at given position."),
          LADA_DECLARE(clear, atomic_center_clear, NOARGS, "Removes all atoms from atomic_center."),
          LADA_DECLARE( extend, atomic_center_extend, O, 
                        "Appends list of atoms to atomic_center.\n\n"
                        "The argument is any iterable objects containing only atoms, "
                        "e.g. another AtomicCenter." ),
          LADA_DECLARE(append, atomic_center_append, O, "Appends an Atom or subtype to the atomic_center.\n"),
          LADA_DECLARE( transform, atomic_center_transform, O, 
                        "Transform a atomic_center in-place.\n\n"
                        "Applies an affine transformation to a atomic_center. "
                        "An affine transformation is a 4x3 matrix, where the upper 3 rows "
                        "correspond to a rotation (applied first), and the last row to "
                        "a translation (applied second).\n\n"
                        ":param matrix:\n"
                        "   Affine transformation defined as a 4x3 numpy array.\n\n"
                        ":returns: A new :class:`AtomicCenter` (or derived) instance." ),
          LADA_DECLARE( __getitem__, atomic_center_subscript, O|METH_COEXIST,
                        "Retrieves atom or slice.\n\n"
                        ":param index:\n"
                        "    If an integer, returns a refence to that atom. "
                            "If a slice, returns a list with all atoms in that slice." ),
          LADA_DECLARE( __setitem__, atomic_center_setitemnormal, VARARGS|METH_COEXIST,
                        "Sets atom or atoms.\n\n"
                        ":param index:\n"
                        "    If an integers, sets that atom to the input value. "
                        "    If a slice, then sets all atoms in refered to in the atomic_center "
                             "to the corresponding atom in the value.\n"
                        ":param value:\n"
                        "    If index is an integer, this should be an :class:`atom <Atom>`. "
                            "If index is a slice, then this should be a sequence of "
                            ":class:`atoms <Atom>` of the exact length of the "
                            "slice."),
          {NULL}  /* Sentinel */
      };
#     undef LADA_DECLARE
  
      static PyTypeObject dummy = {
          PyObject_HEAD_INIT(NULL)
          0,                                 /*ob_size*/
          "lada.vff.cppwrappers.AtomicCenter",   /*tp_name*/
          sizeof(AtomicCenterData),             /*tp_basicsize*/
          0,                                 /*tp_itemsize*/
          (destructor)atomic_center_dealloc,     /*tp_dealloc*/
          0,                                 /*tp_print*/
          0,                                 /*tp_getattr*/
          0,                                 /*tp_setattr*/
          0,                                 /*tp_compare*/
          (reprfunc)atomic_center_repr,          /*tp_repr*/
          0,                                 /*tp_as_number*/
          &atomic_center_as_sequence,            /*tp_as_sequence*/
          &atomic_center_as_mapping,             /*tp_as_mapping*/
          0,                                 /*tp_hash */
          0,                                 /*tp_call*/
          0,                                 /*tp_str*/
          0,                                 /*tp_getattro*/
          0,                                 /*tp_setattro*/
          0,                                 /*tp_as_buffer*/
          Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC | Py_TPFLAGS_HAVE_ITER, /*tp_flags*/
          "Defines a atomic_center.\n\n"         /*tp_doc*/
            "A atomic_center is a special kind of sequence containing only "
            ":class:`Atom`. It also sports attributes such a "
            "cell and scale.\n\n"
            "__init__ accepts different kind of input.\n"
            "  - if 9 numbers are given as arguments, these create the cell vectors, "
                 "where the first three numbers define the first row of the matrix"
                 "(not the first cell vector, eg column).\n"
            "  - if only one object is given it should be the cell matrix.\n"
            "  - the cell can also be given as a keyword argument.\n"
            "  - the scale can only be given as a keyword argument.\n"
            "  - All other keyword arguments become attributes. "
                 "In other words, one could add ``magnetic=0.5`` if one wanted to "
                 "specify the magnetic moment of a atomic_center. It would later be "
                 "accessible as an attribute, eg as ``atomic_center.magnetic``.\n\n"
            ".. note:: The cell is always owned by the object. "
            "Two atomic_centers will not own the same cell object. "
            "The cell given on input is *copied*, *not* referenced. "
            "All other attributes behave like other python attributes: "
            "they are refence if complex objects and copies if a basic python type.",
          (traverseproc)atomic_center_traverse,  /* tp_traverse */
          (inquiry)atomic_center_gcclear,        /* tp_clear */
          0,		                     /* tp_richcompare */
          offsetof(AtomicCenterData, weakreflist),   /* tp_weaklistoffset */
          (getiterfunc)atomic_centeriterator_create,  /* tp_iter */
          0,		                     /* tp_iternext */
          methods,                           /* tp_methods */
          members,                           /* tp_members */
          getsetters,                        /* tp_getset */
          0,                                 /* tp_base */
          0,                                 /* tp_dict */
          0,                                 /* tp_descr_get */
          0,                                 /* tp_descr_set */
          offsetof(AtomicCenterData, pydict),   /* tp_dictoffset */
          (initproc)atomic_center_init,          /* tp_init */
          0,                                 /* tp_alloc */
          (newfunc)PyAtomicCenter_NewWithArgs,  /* tp_new */
      };
      if(dummy.tp_getattro == 0) dummy.tp_getattro = PyObject_GenericGetAttr;
      if(dummy.tp_setattro == 0) dummy.tp_setattro = PyObject_GenericSetAttr;
      return &dummy;
    }

    // Returns pointer to atomic_center iterator type.
    PyTypeObject* atomic_centeriterator_type()
    { 
      static PyTypeObject type = {
          PyObject_HEAD_INIT(NULL)
          0,                                          /*ob_size*/
          "lada.vff.cppwrappers.AtomicCenterIter",   /*tp_name*/
          sizeof(AtomicCenterIterator),                  /*tp_basicsize*/
          0,                                          /*tp_itemsize*/
          (destructor)atomic_centeriterator_dealloc,      /*tp_dealloc*/
          0,                                          /*tp_print*/
          0,                                          /*tp_getattr*/
          0,                                          /*tp_setattr*/
          0,                                          /*tp_compare*/
          0,                                          /*tp_repr*/
          0,                                          
          0,                                          /*tp_as_sequence*/
          0,                                          /*tp_as_mapping*/
          0,                                          /*tp_hash */
          0,                                          /*tp_call*/
          0,                                          /*tp_str*/
          0,                                          /*tp_getattro*/
          0,                                          /*tp_setattro*/
          0,                                          /*tp_as_buffer*/
          Py_TPFLAGS_HAVE_ITER,
          "Iterator over atoms in a atomic_center.",
          0,                                          /* tp_traverse */
          0,                                          /* tp_clear */
          0,                                          /* tp_richcompare */
          0,		                              /* tp_weaklistoffset */
          (getiterfunc)atomic_centeriterator_iter,        /* tp_iter */
          (iternextfunc)atomic_centeriterator_next,       /* tp_iternext */
      };
      return &type;
    }
  } // namespace Crystal
} // namespace LaDa
