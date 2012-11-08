#ifndef LADA_CRYSTAL_HFDATA_H
#define LADA_CRYSTAL_HFDATA_H

#include "LaDaConfig.h"

#include <vector>
#include <ostream>

//! \def PyHFTransform_Check(object)
//!      Returns true if an object is a struture or subtype.
#define PyHFTransform_Check(object) PyObject_TypeCheck(object, LaDa::crystal::hf_type())
//! \def PyHFTransform_CheckExact(object)
//!      Returns true if an object is a hftransform.
#define PyHFTransform_CheckExact(object) object->ob_type == LaDa::crystal::hf_type()
      

namespace LaDa 
{
  namespace crystal
  {
    extern "C" 
    {
      //! \brief Holds data of a hf transform
      //! \details Instances of this object are exactly those that are seen
      //!          within the python interface. C++, however, defines a
      //!          secondary HFTransfirm object which wrapps around a python
      //!          refence to instances of this object. HFTransfrom provides some
      //!          syntactic sugar for handling in c++. 
      struct HFTransformData
      {
        PyObject_HEAD 
        //! The transform to go the hf normal form.
        math::rMatrix3d transform;
        //! Vector of atom wrappers.
        math::iVector3d quotient;
      };
      //! Creates a new hftransform with a given type.
      HFTransformData* PyHFTransform_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      //! Creates a new hftransform with a given type, also calling initialization.
      HFTransformData* PyHFTransform_NewFromArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      //! Creates a deepcopy of hftransform.
      HFTransformData *PyHFTransform_Copy(HFTransformData* _self, PyObject *_memo = NULL);
      // Returns pointer to hftransform type.
      PyTypeObject* hftransform_type();
      //! \brief Initializes a new hftransform from input lattice unit-cell and supercell.
      //! \details Performs initialization from c++ arguments.
      bool hf_transform_init( HFTransformData* _self, 
                                 math::rMatrix3d const &_lattice,
                                 math::rMatrix3d const &_supercell );
    }

  } // namespace Crystal

} // namespace LaDa

#endif
