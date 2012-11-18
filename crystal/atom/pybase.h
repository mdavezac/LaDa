#ifndef LADA_CRYSTAL_ATOMOBJECT_H
#define LADA_CRYSTAL_ATOMOBJECT_H

#include "LaDaConfig.h"

#include <Python.h>

#include <math/eigen.h>

namespace LaDa
{
  namespace crystal
  {

    //! \brief Describes basic atom type. 
    //! \details This is a python object. 
    struct PyAtomObject
    {
      PyObject_HEAD
      //! Holds list of weak pointers.
      PyObject *weakreflist;
      //! Holds python attribute dictionary.
      PyObject *pydict;
      //! Holds the occupation object.
      PyObject *type;
      //! The atomic position in cartesian coordinate.
      math::rVector3d pos;
    };
    
#   ifdef LADA_CRYSTAL_MODULE
      //! Creates a new atom.
      static PyAtomObject* new_atom();
      //! Creates a deepcopy of atom.
      static PyAtomObject *copy_atom(PyAtomObject* _self, PyObject *_memo=NULL);
      //! Creates a new atom with a given type, also calling initialization.
      static PyAtomObject* new_atom(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      //! Returns pointer to atom type.
      static PyTypeObject* atom_type();
      //! Checks type of an object.
      static inline bool check_atom(PyObject *_self)
        { return PyObject_TypeCheck(_self, atom_type()); }
      //! Checks type of an object.
      static inline bool checkexact_atom(PyObject *_self)
        { return _self->ob_type ==  atom_type(); }
      //! Checks type of an object.
      static inline bool check_atom(PyAtomObject *_self)
        { return PyObject_TypeCheck(_self, atom_type()); }
      //! Checks type of an object.
      static inline bool checkexact_atom(PyAtomObject *_self)
        { return _self->ob_type ==  atom_type(); }
#   else
      //! Returns pointer to atom type.
      inline PyTypeObject* atom_type()
        { return (*(PyTypeObject*(*)())api_capsule[0])(); }
      //! Creates a new atom.
      inline PyAtomObject* new_atom()
        { return (*(PyAtomObject*(*)())api_capsule[1])(); }
      //! Creates a new atom with a given type, also calling initialization.
      inline PyAtomObject* new_atom(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs)
        { return (*(PyAtomObject*(*)(PyTypeObject*, PyObject*, PyObject*))
                   api_capsule[2])(_type, _args, _kwargs); }
      //! Creates a deepcopy of atom.
      inline PyAtomObject *copy_atom(PyAtomObject* _self, PyObject *_memo = NULL)
        { return (*(PyAtomObject*(*)(PyAtomObject*, PyObject*))api_capsule[3])(_self, _memo); }
      //! Checks type of an object.
      inline bool check_atom(PyObject *_self)
        { return PyObject_TypeCheck(_self, atom_type()); }
      //! Checks type of an object.
      inline bool check_atom(PyAtomObject *_self)
        { return PyObject_TypeCheck(_self, atom_type()); }
      //! Checks type of an object.
      inline bool checkexact_atom(PyObject *_self)
        { return _self->ob_type ==  atom_type(); }
      //! Checks type of an object.
      inline bool checkexact_atom(PyAtomObject *_self)
        { return _self->ob_type ==  atom_type(); }
#   endif
  }
}
  
#endif
