#include "LaDaConfig.h"

#include <Python.h>
#include <structmember.h>
#define PY_ARRAY_UNIQUE_SYMBOL enumeration_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>

#include <vector>

//! \def PyFCIterator_Check(object)
//!      Returns true if an object is a struture or subtype.
#define PyFCIterator_Check(object) PyObject_TypeCheck(object, LaDa::crystal::structure_type())
//! \def PyFCIterator_CheckExact(object)
//!      Returns true if an object is a structure.
#define PyFCIterator_CheckExact(object) object->ob_type == LaDa::crystal::structure_type()
      

namespace LaDa
{
  namespace enumeration
  {
    //! Type used internally by the counter.
    typedef npy_bool t_fc;
    extern "C" 
    {
      //! \brief Describes basic structure type. 
      //! \details Instances of this object are exactly those that are seen
      //!          within the python interface. C++, however, defines a
      //!          secondary FCIterator object which wrapps around a python
      //!          refence to instances of this object. FCIterator provides some
      //!          syntactic sugar for handling in c++. 
      struct FCIterator
      {
        PyObject_HEAD 
        //! Read-only numpy array referencing inner counter.
        PyArrayObject *yielded;
        //! Inner counter;
        std::vector<t_fc> counter;
        //! Whether this is the first iteration.
        bool is_first;
        //! Concentration
        int ntrue;
      };
      //! Creates a new structure.
      FCIterator* PyFCIterator_New();
      //! Creates a new structure with a given type.
      FCIterator* PyFCIterator_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      //! Creates a new structure with a given type, also calling initialization.
      FCIterator* PyFCIterator_NewFromArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      // Returns pointer to structure type.
      PyTypeObject* fciterator_type();
    }
  }
}
