#if LADA_CRYSTAL_MODULE != 1
#include "LaDaConfig.h"

#include <vector>
#include <ostream>

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
#endif 

#if LADA_CRYSTAL_MODULE != 1
  // Returns pointer to structure type.
  PyTypeObject* structure_type()
    LADA_END( { return (PyTypeObject*)api_capsule[BOOST_PP_SLOT(1)]; } )
#else
  api_capsule[BOOST_PP_SLOT(1)] = (void *)structure_type();
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(BOOST_PP_SLOT(1))
#include BOOST_PP_ASSIGN_SLOT(1)

#if LADA_CRYSTAL_MODULE != 1
  //! Returns address of structure iterator type object.
  PyTypeObject* structureiterator_type()
   LADA_END({ return (PyTypeObject*)api_capsule[BOOST_PP_SLOT(1)]; })
#else
  api_capsule[BOOST_PP_SLOT(1)] = (void *)structureiterator_type();
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(BOOST_PP_SLOT(1))
#include BOOST_PP_ASSIGN_SLOT(1)

#if LADA_CRYSTAL_MODULE != 1
  //! Creates a new structure.
  PyStructureObject* new_structure()
    LADA_END({ return (PyStructureObject*)api_capsule[BOOST_PP_SLOT(1)]; })
#else
  api_capsule[BOOST_PP_SLOT(1)] = (void *)((PyStructureObject*(*)())new_structure);
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(BOOST_PP_SLOT(1))
#include BOOST_PP_ASSIGN_SLOT(1)

#if LADA_CRYSTAL_MODULE != 1
  //! Creates a new structure with a given type, also calling initialization.
  PyStructureObject* new_structure(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs)
    LADA_END( { return (*(PyStructureObject*(*)(PyTypeObject*, PyObject*, PyObject*))
                       api_capsule[BOOST_PP_SLOT(1)])(_type, _args, _kwargs); } )
#else
  api_capsule[BOOST_PP_SLOT(1)] = (void *)((PyStructureObject*(*)(PyTypeObject*, PyObject*, PyObject*))new_structure);
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(BOOST_PP_SLOT(1))
#include BOOST_PP_ASSIGN_SLOT(1)

#if LADA_CRYSTAL_MODULE != 1
  //! Creates a deepcopy of structure.
  PyStructureObject *copy_structure(PyStructureObject* _self, PyObject *_memo=NULL)
    LADA_END( { return (*(PyStructureObject*(*)(PyStructureObject*, PyObject*))
                       api_capsule[BOOST_PP_SLOT(1)])(_self, _memo); } )
#else
  api_capsule[BOOST_PP_SLOT(1)] = (void *)copy_structure;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(BOOST_PP_SLOT(1))
#include BOOST_PP_ASSIGN_SLOT(1)

#if LADA_CRYSTAL_MODULE != 1
  //! Transforms a structure in-place, according to symop.
  void itransform_structure( PyStructureObject* _self,
                             Eigen::Matrix<types::t_real, 4, 3> const &_op )
    LADA_END( { return (*(void(*)(PyStructureObject*, Eigen::Matrix<types::t_real, 4, 3> const &))
                        api_capsule[BOOST_PP_SLOT(1)])(_self, _op); } )
#else
  api_capsule[BOOST_PP_SLOT(1)] = (void *)itransform_structure;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(BOOST_PP_SLOT(1))
#include BOOST_PP_ASSIGN_SLOT(1)

#if LADA_CRYSTAL_MODULE != 1
        //! Checks type of an object.
        bool check_structure(PyObject *_self)
          { return PyObject_TypeCheck(_self, structure_type()); }
        //! Checks type of an object.
        bool checkexact_structure(PyObject *_self)
          { return _self->ob_type ==  structure_type(); }
        //! Checks type of an object.
        bool check_structure(PyStructureObject *_self)
          { return PyObject_TypeCheck(_self, structure_type()); }
        //! Checks type of an object.
        bool checkexact_structure(PyStructureObject *_self)
          { return _self->ob_type ==  structure_type(); }
    } // unamed namespace -- all declarations are static.
  } // namespace Crystal
} // namespace LaDa
#endif
