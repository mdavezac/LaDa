#ifndef LADA_CRYSTAL_ATOM_BASE_H
#define LADA_CRYSTAL_ATOM_BASE_H

#include "LaDaConfig.h"

#define PY_ARRAY_UNIQUE_SYMBOL crystal_ARRAY_API
#define NO_IMPORT_ARRAY

#include <Python.h>

#include <math/eigen.h>

#define LADA_ACQUIRE_PYOBJECT(a, b, type) \
      {                                   \
        PyObject *dummy = (PyObject*)a;   \
        a = (type*)b;                     \
        Py_XINCREF(b);                    \
        Py_XDECREF(dummy);                \
      }

//! Returns true if an object is a string atom or subtype.
#define PyAtom_Check(object) PyObject_TypeCheck(object, LaDa::crystal::atom_type())
//! Returns true if an object is a string atom.
#define PyAtom_CheckExact(object) object->ob_type == LaDa::crystal::atom_type()

namespace LaDa
{
  namespace crystal
  {
    extern "C"
    {
      //! Describes basic atom type. Should be a POD.
      struct AtomData
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
      AtomData* PyAtom_New();
      //! Creates a new atom with a given type.
      AtomData* PyAtom_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      //! Creates a deepcopy of atom.
      AtomData *PyAtom_Copy(AtomData* _self, PyObject *_memo = NULL);
      // Returns pointer to atom type.
      PyTypeObject* atom_type();
    } // extern "C"

    //! \brief imports crystal python module.
    //! \details Sets python exception on import failure.
    bool import();
  } // namespace Crystal
} // namespace LaDa
  
#endif
