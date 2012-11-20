#if LADA_CRYSTAL_MODULE != 1
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
#endif

#if LADA_CRYSTAL_MODULE != 1
  // Returns pointer to hftransform type.
  PyTypeObject* hftransform_type()
    LADA_END({ return (PyTypeObject*)api_capsule[BOOST_PP_SLOT(1)]; })
#else
  api_capsule[BOOST_PP_SLOT(1)] = (void *)hftransform_type();
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(BOOST_PP_SLOT(1))
#include BOOST_PP_ASSIGN_SLOT(1)

#if LADA_CRYSTAL_MODULE != 1
  //! Creates a new hftransform with a given type, also calling initialization.
  PyHFTObject* new_hftransform(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs)
    LADA_END( { return (*(PyHFTObject*(*)(PyTypeObject*, PyObject*, PyObject*))
                        api_capsule[BOOST_PP_SLOT(1)])(_type, _args, _kwargs); } )
#else
  api_capsule[BOOST_PP_SLOT(1)] = (void *)new_hftransform;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(BOOST_PP_SLOT(1))
#include BOOST_PP_ASSIGN_SLOT(1)

#if LADA_CRYSTAL_MODULE != 1
  //! Creates a deepcopy of hftransform.
  PyHFTObject *copy_hftransform(PyHFTObject* _self, PyObject *_memo = NULL)
    LADA_END( { return (*(PyHFTObject*(*)(PyHFTObject*, PyObject*))
                        api_capsule[BOOST_PP_SLOT(1)])(_self, _memo); } )
#else
  api_capsule[BOOST_PP_SLOT(1)] = (void *)copy_hftransform;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(BOOST_PP_SLOT(1))
#include BOOST_PP_ASSIGN_SLOT(1)

#if LADA_CRYSTAL_MODULE != 1
  //! \brief Initializes a new hftransform from input lattice unit-cell and supercell.
  //! \details Performs initialization from c++ arguments.
  bool _init_hft( PyHFTObject* _self, 
                  math::rMatrix3d const &_lattice,
                  math::rMatrix3d const &_supercell )
    LADA_END( { return (*(bool(*)( PyHFTObject*, math::rMatrix3d const&, 
                                   math::rMatrix3d const &))
                        api_capsule[BOOST_PP_SLOT(1)])(_self, _lattice, _supercell); } )
#else
  api_capsule[BOOST_PP_SLOT(1)] = (void *)_init_hft;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(BOOST_PP_SLOT(1))
#include BOOST_PP_ASSIGN_SLOT(1)

#if LADA_CRYSTAL_MODULE != 1
    } // anonymous namespace
  } // namespace Crystal
} // namespace LaDa
#endif
