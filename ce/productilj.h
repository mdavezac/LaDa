#ifndef LADA_ENUM_NDIMITERATOR_H
#define LADA_ENUM_NDIMITERATOR_H
#include "LaDaConfig.h"

#include <Python.h>

#include <vector>

//! \def PyNDimIterator_Check(object)
//!      Returns true if an object is a struture or subtype.
#define PyNDimIterator_Check(object) PyObject_TypeCheck(object, LaDa::crystal::structure_type())
//! \def PyNDimIterator_CheckExact(object)
//!      Returns true if an object is a structure.
#define PyNDimIterator_CheckExact(object) object->ob_type == LaDa::crystal::structure_type()
      

namespace LaDa
{
  namespace ce
  {
    extern "C" 
    {
      //! Product iterator over a sequence, with i < j < ...
      struct ProductILJIterator
      {
        PyObject_HEAD 
        //! Reference to sequence. 
        PyObject *sequence;
        //! Size of the sequence.
        Py_ssize_t N;
        //! Inner counter;
        std::vector<Py_ssize_t> counter;
        //! Whether this is the first iteration.
        bool is_first;
      };
      //! Creates a new structure.
      ProductILJIterator* PyNDimIterator_New();
      //! Creates a new structure with a given type.
      ProductILJIterator* PyNDimIterator_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      //! Creates a new structure with a given type, also calling initialization.
      ProductILJIterator* PyNDimIterator_NewFromArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      // Returns pointer to structure type.
      PyTypeObject* productiljiterator_type();
    }
  }
}
#endif
