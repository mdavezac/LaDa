#ifndef LADA_CRYSTAL_ATOMOBJECT_H
#define LADA_CRYSTAL_ATOMOBJECT_H

#include "LaDaConfig.h"

#include <Python.h>

#include <math/eigen.h>

namespace LaDa
{
  namespace crystal
  {
    namespace // unamed namespace -- declarations unavailable outside of scope.
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
      
#     ifdef LADA_CRYSTAL_MODULE
        //! Returns pointer to atom type.
        PyTypeObject* atom_type();
#       define LADA_VALUE ((PyAtomObject*(*)())new_atom)
#       include "../_incrementor.hpp"
        //! Creates a new atom.
        PyAtomObject* new_atom();
#       define LADA_VALUE ((PyAtomObject*(*)())new_atom)
#       include "../_incrementor.hpp"
        //! Creates a new atom with a given type, also calling initialization.
        PyAtomObject* new_atom(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
#       define LADA_VALUE ((PyAtomObject*(*)(PyTypeObject*, PyObject*, PyObject*))new_atom)
#       include "../_incrementor.hpp"
        //! Creates a deepcopy of atom.
        PyAtomObject *copy_atom(PyAtomObject* _self, PyObject *_memo=NULL);
#       define LADA_VALUE copy_atom 
#       include "../_incrementor.hpp"
        //! Checks type of an object.
        inline bool check_atom(PyObject *_self)
          { return PyObject_TypeCheck(_self, atom_type()); }
        //! Checks type of an object.
        inline bool checkexact_atom(PyObject *_self)
          { return _self->ob_type ==  atom_type(); }
        //! Checks type of an object.
        inline bool check_atom(PyAtomObject *_self)
          { return PyObject_TypeCheck(_self, atom_type()); }
        //! Checks type of an object.
        inline bool checkexact_atom(PyAtomObject *_self)
          { return _self->ob_type ==  atom_type(); }
#     else
        //! Returns pointer to atom type.
        inline PyTypeObject* atom_type()
          { return (PyTypeObject*)api_capsule[BOOST_PP_SLOT(1)]; }
#       define LADA_VALUE atom_type()
#       include "../_incrementor.hpp"
        //! Creates a new atom.
        inline PyAtomObject* new_atom()
          { return (*(PyAtomObject*(*)())api_capsule[BOOST_PP_SLOT(1)])(); }
#       define LADA_VALUE ((PyAtomObject*(*)())new_atom)
#       include "../_incrementor.hpp"
        //! Creates a new atom with a given type, also calling initialization.
        inline PyAtomObject* new_atom(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs)
          { return (*(PyAtomObject*(*)(PyTypeObject*, PyObject*, PyObject*))
                     api_capsule[BOOST_PP_SLOT(1)])(_type, _args, _kwargs); }
#       define LADA_VALUE ((PyAtomObject*(*)(PyTypeObject*, PyObject*, PyObject*))new_atom)
#       include "../_incrementor.hpp"
        //! Creates a deepcopy of atom.
        inline PyAtomObject *copy_atom(PyAtomObject* _self, PyObject *_memo = NULL)
          { return (*(PyAtomObject*(*)(PyAtomObject*, PyObject*))api_capsule[BOOST_PP_SLOT(1)])(_self, _memo); }
#       define LADA_VALUE copy_atom 
#       include "../_incrementor.hpp"
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
#     endif
    } // static namespace
  }
}
  
#endif
