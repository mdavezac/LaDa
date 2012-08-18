#ifndef LADA_ENUM_NDIMITERATOR_H
#define LADA_ENUM_NDIMITERATOR_H
#include "LaDaConfig.h"

#include <Python.h>
#include <structmember.h>
#define PY_ARRAY_UNIQUE_SYMBOL enumeration_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>

#include <vector>

//! \def PyNDimIterator_Check(object)
//!      Returns true if an object is a struture or subtype.
#define PyNDimIterator_Check(object) PyObject_TypeCheck(object, LaDa::crystal::structure_type())
//! \def PyNDimIterator_CheckExact(object)
//!      Returns true if an object is a structure.
#define PyNDimIterator_CheckExact(object) object->ob_type == LaDa::crystal::structure_type()
      

namespace LaDa
{
  namespace enumeration
  {
    //! Type used internally by the counter.
    typedef npy_short t_ndim;
    extern "C" 
    {
      //! \brief Describes basic structure type. 
      //! \details Instances of this object are exactly those that are seen
      //!          within the python interface. C++, however, defines a
      //!          secondary NDimIterator object which wrapps around a python
      //!          refence to instances of this object. NDimIterator provides some
      //!          syntactic sugar for handling in c++. 
      struct NDimIterator
      {
        PyObject_HEAD 
        //! Holds beginning and end of range.
        std::vector<t_ndim> ends;
        //! Read-only numpy array referencing inner counter.
        PyArrayObject *yielded;
        //! Inner counter;
        std::vector<t_ndim> counter;
      };
      //! Creates a new structure.
      NDimIterator* PyNDimIterator_New();
      //! Creates a new structure with a given type.
      NDimIterator* PyNDimIterator_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      //! Creates a new structure with a given type, also calling initialization.
      NDimIterator* PyNDimIterator_NewFromArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      // Returns pointer to structure type.
      PyTypeObject* ndimiterator_type();
    }
  }
}
#endif
