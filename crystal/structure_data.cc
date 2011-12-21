#include "LaDaConfig.h"

#include <Python.h>
#include "structmember.h"

#define PY_ARRAY_UNIQUE_SYMBOL crystal_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>


#include <math/python/python.hpp>
#include <python/exceptions.h>
#include <python/numpy_types.h>

#include "structure_data.h"

// structure getters/settters.
#include "structure_getset.hpp"
// structure member functions.
#include "structure_members.hpp"
// creation, deallocation, initialization.
#include "structure_cdi.hpp"

namespace LaDa
{
  namespace crystal
  {
    //! returns 1 angstrom.
    static PyObject* get_unit_angstrom();
    //! Creates a new structure.
    StructureData* PyStructure_New()
    {
      StructureData* result = (StructureData*) structure_type()->tp_alloc(structure_type(), 0);
      if(not result) return NULL;
      result->weakreflist = NULL;
      result->scale = 1e0;
      new(&result->cell) LaDa::math::rMatrix3d(LaDa::math::rMatrix3d::Identity());
      result->pydict = PyDict_New();
      if(result->pydict == NULL) { Py_DECREF(result); return NULL; }
      return result;
    }
    //! Creates a new structure with a given type.
    StructureData* PyStructure_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs)
    {
      StructureData* result = (StructureData*)_type->tp_alloc(_type, 0);
      if(not result) return NULL;
      result->weakreflist = NULL;
      result->scale = 1e0;
      new(&result->cell) LaDa::math::rMatrix3d(LaDa::math::rMatrix3d::Identity());
      result->pydict = PyDict_New();
      if(result->pydict == NULL) { Py_DECREF(result); return NULL; }
      return result;
    }

    // Creates a deepcopy of structure.
    StructureData *PyStructure_Copy(StructureData* _self, PyObject *_memo)
    {
      StructureData* result = (StructureData*)_self->ob_type->tp_alloc(_self->ob_type, 0);
      if(not result) return NULL;
      result->weakreflist = NULL;
      new(&result->cell) LaDa::math::rMatrix3d(_self->cell);
      result->pydict = NULL;
      result->scale = _self->scale;
      PyObject* copymod = PyImport_ImportModule("copy");
      if(copymod == NULL) return NULL;
      PyObject *deepcopystr = PyString_FromString("deepcopy");
      if(not deepcopystr) { Py_DECREF(copymod); return NULL; }
      else if(_self->pydict != NULL)
      {
        if(_memo == NULL)
          result->pydict = PyObject_CallMethodObjArgs(copymod, deepcopystr, _self->pydict, NULL);
        else
          result->pydict = PyObject_CallMethodObjArgs(copymod, deepcopystr, _self->pydict, _memo, NULL);
        if(result->pydict == NULL) { Py_DECREF(result);  result == NULL; }
      }

      Py_DECREF(copymod);
      Py_DECREF(deepcopystr);
      return result;
    }
    // Returns pointer to structure type.
    PyTypeObject* structure_type()
    {
#     ifdef LADA_DECLARE
#       error LADA_DECLARE already defined.
#     endif
#     define LADA_DECLARE(name, doc) \
        { const_cast<char*>(#name), (getter) structure_get ## name, \
          (setter) structure_set ## name, const_cast<char*>(doc) }
      
      static PyGetSetDef getsetters[] = {
          LADA_DECLARE(cell, "Cell matrix in cartesian coordinates.\n\n"
                             "Unlike most ab-initio codes, cell-vectors are "
                             "given in column vector format. " 
                             "The cell does not yet have units. "
                             "Units depend upon `lada.crystal.Structure.scale`. "
                             "Across lada, it is expected that a cell time "
                             "this scale are angstroms. Finally, the cell "
                             "is owned internally by the structure. It cannot be set to "
                             "reference an object (say a list or numpy array). "
                             "``structure.cell = some_list`` will copy the values of ``some_list``. "),
          LADA_DECLARE(scale, "Scale factor of this structure.\n"
                              "Should be a number or unit given by "
                              "the python package *quantities*. In that case, "
                              "it is converted to angstroms. "),
          {NULL}  /* Sentinel */
      };
#     undef LADA_DECLARE
#     define LADA_DECLARE(name, object, doc) \
        { const_cast<char*>(#name), T_OBJECT_EX, \
          offsetof(StructureData, object), 0, const_cast<char*>(doc) }
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
          LADA_DECLARE(copy, structure_copy, NOARGS, "Returns a deepcopy of the structure."),
          LADA_DECLARE( to_dict, structure_to_dict, NOARGS, 
                        "Returns a dictionary with shallow copies of items." ),
          LADA_DECLARE(__copy__, structure_shallowcopy, NOARGS, "Shallow copy of an structure."),
          LADA_DECLARE(__deepcopy__, PyStructure_Copy, O, "Deep copy of an structure."),
          LADA_DECLARE(__getstate__, structure_getstate, NOARGS, "Implements pickle protocol."),
          LADA_DECLARE(__setstate__, structure_setstate, O, "Implements pickle protocol."),
          LADA_DECLARE(__reduce__,   structure_reduce, NOARGS, "Implements pickle protocol."),
          LADA_DECLARE( add_atom,   structure_add_atom, KEYWORDS,
                        "Adds atom to structure.\n\n"
                        "The argument to this function is either another atom, "
                        "in which case a reference to that atom is appended to "
                        "the structure. Or, it is any arguments used to "
                        "initialize atoms in`lada.crystal.cppwrappers.Atom.__init__`."
                        "Finally, this function can be chained as follows:\n\n"
                        ">>> structure.add_atom(0,0,0, 'Au')\\\n"
                        "...                   (0.25, 0.25, 0.25, ['Pd', 'Si'], m=5)\\\n"
                        "...                   (atom_from_another_structure)\n\n"
                        "In the example above, both ``structure`` and the *other* structure will "
                        "reference the same atom (``atom_from_another_structure``). "
                        "Changing, say, that atom's type in one structure will also "
                        "change it in the other."),
          {NULL}  /* Sentinel */
      };
#     undef LADA_DECLARE
  
      static PyTypeObject dummy = {
          PyObject_HEAD_INIT(NULL)
          0,                                 /*ob_size*/
          "lada.crystal.cppwrappers.Structure",   /*tp_name*/
          sizeof(StructureData),             /*tp_basicsize*/
          0,                                 /*tp_itemsize*/
          (destructor)structure_dealloc,     /*tp_dealloc*/
          0,                                 /*tp_print*/
          0,                                 /*tp_getattr*/
          0,                                 /*tp_setattr*/
          0,                                 /*tp_compare*/
          (reprfunc)structure_repr,          /*tp_repr*/
          0,                                 /*tp_as_number*/
          0,                                 /*tp_as_sequence*/
          0,                                 /*tp_as_mapping*/
          0,                                 /*tp_hash */
          0,                                 /*tp_call*/
          0,                                 /*tp_str*/
          0,                                 /*tp_getattro*/
          0,                                 /*tp_setattro*/
          0,                                 /*tp_as_buffer*/
          Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC, /*tp_flags*/
          "Defines a structure.\n\n"
            "__init__ accepts different kind of input.\n"
            "  - if 9 numbers are given as arguments, these create the cell vectors, "
                 "where the first three numbers define the first row of the matrix"
                 "(not the first cell vector, eg column\n"
            "  - if only one object is given it should be the cell matrix.\n"
            "  - the cell can also be given as a keyword argument.\n"
            "  - the scale can only be given as a keyword argument.\n"
            "  - All other keyword arguments become attributes. "
                 "In other words, one could add ``magnetic=0.5`` if one wanted to "
                 "specify the magnetic moment of a structure.\n\n"
            "Note that the cell is always owned by the object. "
            "Two structures will not own the same cell object. "
            "The cell given on input is *copied*, *not* referenced. "
            "All other attributes behave like other python attributes: "
            "they are refence if complex objects and copies if a basic python type.",
          (traverseproc)structure_traverse,  /* tp_traverse */
          (inquiry)structure_gcclear,        /* tp_clear */
          0,		                     /* tp_richcompare */
          offsetof(StructureData, weakreflist),   /* tp_weaklistoffset */
          0,		                     /* tp_iter */
          0,		                     /* tp_iternext */
          methods,                           /* tp_methods */
          members,                           /* tp_members */
          getsetters,                        /* tp_getset */
          0,                                 /* tp_base */
          0,                                 /* tp_dict */
          0,                                 /* tp_descr_get */
          0,                                 /* tp_descr_set */
          offsetof(StructureData, pydict),   /* tp_dictoffset */
          (initproc)structure_init,          /* tp_init */
          0,                                 /* tp_alloc */
          (newfunc)PyStructure_NewWithArgs,  /* tp_new */
      };
      if(dummy.tp_getattro == 0) dummy.tp_getattro = PyObject_GenericGetAttr;
      if(dummy.tp_setattro == 0) dummy.tp_setattro = PyObject_GenericSetAttr;
      return &dummy;
    }
  } // namespace Crystal
} // namespace LaDa
