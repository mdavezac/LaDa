#ifndef LADA_ENUM_MANIPULATIONS_H
#define LADA_ENUM_MANIPULATIONS_H
#include "LaDaConfig.h"

#include <Python.h>
#include <structmember.h>
#define PY_ARRAY_UNIQUE_SYMBOL enumeration_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>

#include "ndimiterator.h"

//! \def PyManipulations_Check(object)
//!      Returns true if an object is a struture or subtype.
#define PyManipulations_Check(object) PyObject_TypeCheck(object, LaDa::crystal::structure_type())
//! \def PyManipulations_CheckExact(object)
//!      Returns true if an object is a structure.
#define PyManipulations_CheckExact(object) object->ob_type == LaDa::crystal::structure_type()
      

namespace LaDa
{
  namespace enumeration
  {
    //! Type used internally by the counter.
    extern "C" 
    {
      //! \brief Describes basic manipulations iterator. 
      struct Manipulations
      {
        PyObject_HEAD 
        //! Array object
        PyArrayObject *arrayout;
        PyObject *arrayin;
        //! Inner counter;
        std::vector<t_ndim> counter;
        //! Inner counter;
        std::vector<t_ndim> substituted;
        //! Inner counter;
        std::vector< std::vector<t_ndim> > substitutions;
        //! iterator over substitutions
        std::vector< std::vector<t_ndim> >:: iterator i_first;
        //! iterator over substitutions
        std::vector< std::vector<t_ndim> >:: iterator i_end;
        //! Whether this is the first iteration.
        bool is_first;
      };
      //! Creates a new structure.
      Manipulations* PyManipulations_New();
      //! Creates a new structure with a given type.
      Manipulations* PyManipulations_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      //! Creates a new structure with a given type, also calling initialization.
      Manipulations* PyManipulations_NewFromArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      // Returns pointer to structure type.
      PyTypeObject* manipulations_type();
    }
  }
}
#endif
