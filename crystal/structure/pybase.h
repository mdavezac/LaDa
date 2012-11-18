#ifndef LADA_CRYSTAL_PYSTRUCTURE_H
#define LADA_CRYSTAL_PYSTRUCTURE_H

#include "LaDaConfig.h"

#include <vector>
#include <ostream>

#include "../atom/atom.h"

namespace LaDa 
{
  namespace crystal
  {
    namespace // unamed namespace == all declarations are static.
    {
      //! \brief Describes basic structure type. 
      //! \details Instances of this object are exactly those that are seen
      //!          within the python interface. C++, however, defines a
      //!          secondary Structure object which wrapps around a python
      //!          refence to instances of this object. Structure provides some
      //!          syntactic sugar for handling in c++. 
      struct PyStructureObject
      {
        PyObject_HEAD 
        //! Holds list of weak pointers.
        PyObject *weakreflist;
        //! Holds python attribute dictionary.
        PyObject *pydict;
        //! Holds python attribute dictionary.
        PyObject *scale;
        //! The unit-cell of the structure in cartesian coordinate.
        math::rMatrix3d cell;
        //! Vector of atom wrappers.
        std::vector<Atom> atoms;
      };
#     ifdef LADA_CRYSTAL_MODULE
        // Returns pointer to structure type.
        PyTypeObject* structure_type();
        //! Returns address of structure iterator type object.
        PyTypeObject* structureiterator_type();
        //! Creates a new structure.
        PyStructureObject* new_structure();
        //! Creates a new structure with a given type, also calling initialization.
        PyStructureObject* new_structure(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
        //! Creates a deepcopy of structure.
        PyStructureObject *copy_structure(PyStructureObject* _self, PyObject *_memo=NULL);
        //! Transforms a structure in-place, according to symop.
        void itransform_structure( PyStructureObject* _self,
                                          Eigen::Matrix<types::t_real, 4, 3> const &_op );
        //! Checks type of an object.
        inline bool check_structure(PyObject *_self)
          { return PyObject_TypeCheck(_self, structure_type()); }
        //! Checks type of an object.
        inline bool checkexact_structure(PyObject *_self)
          { return _self->ob_type ==  structure_type(); }
        //! Checks type of an object.
        inline bool check_structure(PyStructureObject *_self)
          { return PyObject_TypeCheck(_self, structure_type()); }
        //! Checks type of an object.
        inline bool checkexact_structure(PyStructureObject *_self)
          { return _self->ob_type ==  structure_type(); }
#     else
        //! Returns pointer to structure type.
        inline PyTypeObject* structure_type()
          { return (PyTypeObject*)api_capsule[4]; }
        //! Returns pointer to structure type.
        inline PyTypeObject* structureiterator_type()
          { return (PyTypeObject*)api_capsule[5]; }
        //! Creates a new structure.
        inline PyStructureObject* new_structure()
          { return (*(PyStructureObject*(*)())api_capsule[6])(); }
        //! Creates a new structure with a given type, also calling initialization.
        inline PyStructureObject* new_structure(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs)
          { return (*(PyStructureObject*(*)(PyTypeObject*, PyObject*, PyObject*))
                     api_capsule[7])(_type, _args, _kwargs); }
        //! Transforms a structure in-place, according to symop.
        inline void itransform_structure( PyStructureObject* _self,
                                          Eigen::Matrix<types::t_real, 4, 3> const &_op )
          { return (*(void(*)(PyStructureObject*, Eigen::Matrix<types::t_real, 4, 3> const &))
                     api_capsule[8])(_self, _op); }
        //! Creates a deepcopy of structure.
        inline PyStructureObject *copy_structure(PyStructureObject* _self, PyObject *_memo = NULL)
          { return (*(PyStructureObject*(*)(PyStructureObject*, PyObject*))api_capsule[8])(_self, _memo); }
        //! Checks type of an object.
        inline bool check_structure(PyObject *_self)
          { return PyObject_TypeCheck(_self, structure_type()); }
        //! Checks type of an object.
        inline bool checkexact_structure(PyObject *_self)
          { return _self->ob_type ==  structure_type(); }
        //! Checks type of an object.
        inline bool check_structure(PyStructureObject *_self)
          { return PyObject_TypeCheck(_self, structure_type()); }
        //! Checks type of an object.
        inline bool checkexact_structure(PyStructureObject *_self)
          { return _self->ob_type ==  structure_type(); }
#     endif
    } // unamed namespace -- all declarations are static.
  } // namespace Crystal
} // namespace LaDa

#endif
