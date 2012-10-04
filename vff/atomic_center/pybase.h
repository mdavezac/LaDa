#ifndef LADA_VFF_ATOMIC_CENTER_H
#define LADA_VFF_ATOMIC_CENTER_H

#include "LaDaConfig.h"

#define PY_ARRAY_UNIQUE_SYMBOL lada_vff_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>

#include <vector>
#include <ostream>

#include "../atom/atom.h"

//! \def PyStructure_Check(object)
//!      Returns true if an object is an atomic_center or subtype.
#define PyAtomicCenter_Check(object) PyObject_TypeCheck(object, LaDa::vff::atomic_center_type())
//! \def PyAtomicCenter_CheckExact(object)
//!      Returns true if an object is a atomic_center.
#define PyAtomicCenter_CheckExact(object) object->ob_type == LaDa::vff::atomic_center_type()
      

namespace LaDa 
{
  namespace vff
  {
    extern "C" 
    {
      //! \brief Describes basic atomic_center type. 
      //! \details An atomic center is a node on the first neighbor graph of a
      //!          zinc-blende structure. It holds a reference to a node atom.
      //!          It also holds references to other atomic center which
      //!          represent the first neighbors.
      struct AtomicCenterData
      {
        PyObject_HEAD 
        //! Holds list of weak pointers.
        PyObject *weakreflist;
        //! Holds reference to other bonds.
        std::vector<AtomicCenterData*> bonds;
        //! Holds reference to other an atom.
        crystal::AtomData* center;
        //! Index of the atom in the structure.
        types::t_unsigned index;
      };
      //! Creates a new atomic_center.
      AtomicCenterData* PyAtomicCenter_New();
      //! Creates a new atomic_center with a given type.
      AtomicCenterData* PyAtomicCenter_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      //! Creates a new atomic_center with a given type, also calling initialization.
      AtomicCenterData* PyAtomicCenter_NewFromArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      //! Creates a deepcopy of atomic_center.
      AtomicCenterData *PyAtomicCenter_Copy(AtomicCenterData* _self, PyObject *_memo = NULL);
      // Returns pointer to atomic_center type.
      PyTypeObject* atomic_center_type();
      //! Returns address of atomic_center iterator type object.
      PyTypeObject* atomic_centeriterator_type();
    }

  } // namespace vff

} // namespace LaDa

#endif
