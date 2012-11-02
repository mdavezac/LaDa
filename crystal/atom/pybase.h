#ifndef LADA_CRYSTAL_ATOM_BASE_H
#define LADA_CRYSTAL_ATOM_BASE_H

#include "LaDaConfig.h"

#include <Python.h>

#include <math/eigen.h>

//! Returns true if an object is an atom or subtype.
#define PyAtom_Check(object) PyObject_TypeCheck(object, LaDa::crystal::atom_type())
//! Returns true if an object is an atom.
#define PyAtom_CheckExact(object) object->ob_type == LaDa::crystal::atom_type()

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
    
    //! Creates a new atom.
    PyAtomObject* PyAtom_New();
    //! Creates a new atom with a given type.
    PyAtomObject* PyAtom_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
    //! Creates a new atom with a given type, also calling initialization.
    PyAtomObject* PyAtom_NewFromArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
    //! Creates a deepcopy of atom.
    PyAtomObject *PyAtom_Copy(PyAtomObject* _self, PyObject *_memo = NULL);
    // Returns pointer to atom type.
    PyTypeObject* atom_type();
  }
}
  
#endif
