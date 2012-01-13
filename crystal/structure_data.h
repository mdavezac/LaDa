#ifndef LADA_CRYSTAL_STRUCTUREDATA_H
#define LADA_CRYSTAL_STRUCTUREDATA_H

#include "LaDaConfig.h"

#include <vector>
#include <ostream>

#include <python/quantity.h>
#include "atom.h"

//! Returns true if an object is a struture or subtype.
#define PyStructure_Check(object) PyObject_TypeCheck(object, LaDa::crystal::structure_type())
//! Returns true if an object is a structure.
#define PyStructure_CheckExact(object) object->ob_type == LaDa::crystal::structure_type()
      

namespace LaDa 
{
  namespace crystal
  {
    extern "C" 
    {
      //! \brief Describes basic structure type. 
      //! \details Instances of this object are exactly those that are seen
      //!          within the python interface. C++, however, defines a
      //!          secondary Structure object which wrapps around a python
      //!          refence to instances of this object. Structure provides some
      //!          syntactic sugar for handling in c++. 
      struct StructureData
      {
        PyObject_HEAD 
        //! Holds list of weak pointers.
        PyObject *weakreflist;
        //! Holds python attribute dictionary.
        PyObject *pydict;
        //! The scale in which cartesian units are given. Units of ansgtroms.
        types::t_real scale;
        //! The unit-cell of the structure in cartesian coordinate.
        math::rMatrix3d cell;
        //! Vector of atom wrappers.
        std::vector<Atom> atoms;
      };
      //! Creates a new structure.
      StructureData* PyStructure_New();
      //! Creates a new structure with a given type.
      StructureData* PyStructure_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      //! Creates a new structure with a given type, also calling initialization.
      StructureData* PyStructure_NewFromArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      //! Creates a deepcopy of structure.
      StructureData *PyStructure_Copy(StructureData* _self, PyObject *_memo = NULL);
      // Returns pointer to structure type.
      PyTypeObject* structure_type();
      //! Returns address of structure iterator type object.
      PyTypeObject* structureiterator_type();
    }

  } // namespace Crystal

} // namespace LaDa

#endif
