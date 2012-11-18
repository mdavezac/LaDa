#ifndef LADA_CRYSTAL_HFDATA_H
#define LADA_CRYSTAL_HFDATA_H

#include "LaDaConfig.h"

#include <vector>
#include <ostream>

namespace LaDa 
{
  namespace crystal
  {
    namespace // limits declarations to current translation units.
    {
      //! \brief Holds data of a hf transform
      //! \details Instances of this object are exactly those that are seen
      //!          within the python interface. C++, however, defines a
      //!          secondary HFTransfirm object which wrapps around a python
      //!          refence to instances of this object. HFTransfrom provides some
      //!          syntactic sugar for handling in c++. 
      struct PyHFTObject
      {
        PyObject_HEAD 
        //! The transform to go the hf normal form.
        math::rMatrix3d transform;
        //! Vector of atom wrappers.
        math::iVector3d quotient;
      };
#     ifdef LADA_CRYSTAL_MODULE
        // Returns pointer to hftransform type.
        PyTypeObject* hftransform_type();
        //! Creates a new hftransform with a given type, also calling initialization.
        PyHFTObject* new_hftransform(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
        //! Creates a deepcopy of hftransform.
        PyHFTObject *copy_hftransform(PyHFTObject* _self, PyObject *_memo = NULL);
        //! \brief Initializes a new hftransform from input lattice unit-cell and supercell.
        //! \details Performs initialization from c++ arguments.
        bool _init_hft( PyHFTObject* _self, 
                        math::rMatrix3d const &_lattice,
                        math::rMatrix3d const &_supercell );
#     else
        // Returns pointer to hftransform type.
        inline PyTypeObject* hftransform_type()
          { return (PyTypeObject*)api_capsule[20]; }
        //! Creates a new hftransform with a given type, also calling initialization.
        inline PyHFTObject* new_hftransform(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs)
          { return (*(PyHFTObject*(*)(PyTypeObject*, PyObject*, PyObject*))
                     api_capsule[21])(_type, _args, _kwargs); }
        //! Creates a deepcopy of hftransform.
        inline PyHFTObject *copy_hftransform(PyHFTObject* _self, PyObject *_memo = NULL)
          { return (*(PyHFTObject*(*)(PyHFTObject*, PyObject*))
                     api_capsule[22])(_self, _memo); }
        //! \brief Initializes a new hftransform from input lattice unit-cell and supercell.
        //! \details Performs initialization from c++ arguments.
        inline bool _init_hft( PyHFTObject* _self, 
                               math::rMatrix3d const &_lattice,
                               math::rMatrix3d const &_supercell )
          { return (*(bool(*)( PyHFTObject*, math::rMatrix3d const&, 
                               math::rMatrix3d const &))
                     api_capsule[23])(_self, _lattice, _supercell); }
#     endif
    } // anonymous namespace
  } // namespace Crystal
} // namespace LaDa

#endif
