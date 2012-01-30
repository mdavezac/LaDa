#ifndef LADA_CRYSTAL_SMITHDATA_H
#define LADA_CRYSTAL_SMITHDATA_H

#include "LaDaConfig.h"

#define PY_ARRAY_UNIQUE_SYMBOL crystal_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>

#include <vector>
#include <ostream>

//! \def PySmithTransform_Check(object)
//!      Returns true if an object is a struture or subtype.
#define PySmithTransform_Check(object) PyObject_TypeCheck(object, LaDa::crystal::smith_type())
//! \def PySmithTransform_CheckExact(object)
//!      Returns true if an object is a smithtransform.
#define PySmithTransform_CheckExact(object) object->ob_type == LaDa::crystal::smith_type()
      

namespace LaDa 
{
  namespace crystal
  {
    extern "C" 
    {
      //! \brief Holds data of a smith transform
      //! \details Instances of this object are exactly those that are seen
      //!          within the python interface. C++, however, defines a
      //!          secondary SmithTransfirm object which wrapps around a python
      //!          refence to instances of this object. SmithTransfrom provides some
      //!          syntactic sugar for handling in c++. 
      struct SmithTransformData
      {
        PyObject_HEAD 
        //! Holds python attribute dictionary.
        PyObject *pydict;
        //! The transform to go the smith normal form.
        math::rMatrix3d transform;
        //! Vector of atom wrappers.
        math::iVector3d quotient;
      };
      //! Creates a new smithtransform with a given type.
      SmithTransformData* PySmithTransform_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      //! Creates a new smithtransform with a given type, also calling initialization.
      SmithTransformData* PySmithTransform_NewFromArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      //! Creates a deepcopy of smithtransform.
      SmithTransformData *PySmithTransform_Copy(SmithTransformData* _self, PyObject *_memo = NULL);
      // Returns pointer to smithtransform type.
      PyTypeObject* smithtransform_type();
      //! \brief Initializes a new smithtransform from input lattice unit-cell and supercell.
      //! \details Performs initialization from c++ arguments.
      bool smith_transform_init( SmithTransformData* _self, 
                                 math::rMatrix3d const &_lattice,
                                 math::rMatrix3d const &_supercell );
    }

  } // namespace Crystal

} // namespace LaDa

#endif
