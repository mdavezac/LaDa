#ifndef LADA_PYTHON_PYTHON_H
#define LADA_PYTHON_PYTHON_H
#ifndef __cplusplus
# error LaDa requires a cpp compiler
#endif

#ifndef LADA_PYTHON_MODULE
#  define LADA_PYTHON_MODULE 100
#endif

#if LADA_PYTHON_MODULE != 1

# include "LaDaConfig.h"
# include <Python.h>
# include <structmember.h>

# include <vector>
# include <string>
# include <iostream>

# include <boost/mpl/int.hpp>
# include <boost/type_traits/is_floating_point.hpp>
# include <boost/preprocessor/arithmetic/inc.hpp>
# include <python/ppslot.hpp>
# define BOOST_PP_VALUE 0
# include LADA_ASSIGN_SLOT(python)

# include <Eigen/Core>

# include <errors/exceptions.h>
# include "types.h"

  namespace LaDa
  {
    namespace python
    {

#     if LADA_PYTHON_MODULE == 100
        /* This section is used in modules that use lada.python's API */
#       ifdef LADA_NO_IMPORT
          extern
#       endif 
        void **api_capsule;
        
        namespace 
        {
          // Return -1 on error, 0 on success.
          // PyCapsule_Import will set an exception if there's an error.
          inline bool import(void)
          {
            PyObject *module = PyImport_ImportModule("lada.cppwrappers");
            if(not module) return false;
            Py_DECREF(module);
            api_capsule = (void **)PyCapsule_Import("lada.cppwrappers._C_API", 0);
            return api_capsule != NULL;
          }
        }
#     endif
#else
# define BOOST_PP_VALUE 0
# include LADA_ASSIGN_SLOT(python)
#endif

#if LADA_PYTHON_MODULE != 1
#  ifdef LADA_INLINE
#    error LADA_INLINE already defined
#  endif
#  if LADA_PYTHON_MODULE == 100
#    define LADA_INLINE inline
#  elif LADA_PYTHON_MODULE == 0
#    define LADA_INLINE
#  endif
#  ifdef LADA_END
#    error LADA_END already defined
#  elif LADA_PYTHON_MODULE == 0
#    define LADA_END(X) ;
#  elif LADA_PYTHON_MODULE == 100
#    define LADA_END(X) { X }
#  endif
#endif

#if LADA_PYTHON_MODULE != 1
  class Object;
  namespace
  {
#endif
#if LADA_PYTHON_MODULE != 1
  //! Object reset function.
  //! Declared as friend to object so that it can be linked at runtime.
  LADA_INLINE void object_reset(PyObject*& _object, PyObject *_in)
    LADA_END( return ( *(void(*)(PyObject*&, PyObject*))
                       api_capsule[LADA_SLOT(python)])(_object, _in); ) 
#else
  api_capsule[LADA_SLOT(python)] = (void *)((void(*)(PyObject*&, PyObject*)) object_reset);
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(python))
#include LADA_ASSIGN_SLOT(python)
  
#if LADA_PYTHON_MODULE != 1
  LADA_INLINE bool object_equality_op(Object const& _self, Object const &_b)
    LADA_END( return ( *(bool(*)(Object const&, Object const&))
                       api_capsule[LADA_SLOT(python)])(_self, _b); )
#else
  api_capsule[LADA_SLOT(python)] = (void *)object_equality_op;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(python))
#include LADA_ASSIGN_SLOT(python)

#if LADA_PYTHON_MODULE != 1
  LADA_INLINE std::ostream& operator<<(std::ostream &_stream, Object const &_ob)
    LADA_END( return ( (std::ostream&(*)(std::ostream&, Object const&)) 
                       api_capsule[LADA_SLOT(python)] )(_stream, _ob); )
#else
  api_capsule[LADA_SLOT(python)] = (void*)
       ( (std::ostream&(*)(std::ostream&, Object const&)) operator<< );
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(python))
#include LADA_ASSIGN_SLOT(python)

#if LADA_PYTHON_MODULE != 1
  }
#endif

// in namespace LaDa::python, but not anonymous.
# include "object.h"

#if LADA_PYTHON_MODULE != 1
  namespace
  {
    //! \brief Acquires a reference to an object.
    //! \details Input is XINCREF'ed before the return wrappers is created. 
    inline Object acquire(PyObject *_in) { Py_XINCREF(_in); return Object(_in); }
    //! \brief Steals a reference to an object.
    //! \details Input is XINCREF'ed before the return wrappers is created. 
    inline Object steal(PyObject *_in) { return Object(_in); }
    
    //! \brief Dumps representation of an object.
    //! \details Will throw c++ exceptions if python calls fail. Does not clear
    //!          python exceptions.
    inline std::ostream& operator<< (std::ostream &stream, PyObject* _ob)
      { return stream << Object::acquire(_ob); }
  }
#endif

// in namespace LaDa::python, but not anonymous.
# include "random_access_list_iterator.h"
# include "random_access_tuple_iterator.h"

#if LADA_PYTHON_MODULE != 1
  namespace numpy
  {
    namespace 
    {
#endif

#     include "numpy_types.h"
#     include "wrap_numpy.h"

#if LADA_PYTHON_MODULE != 1
    }
  }
  namespace 
  {
#endif

# include "quantity.h"

#if LADA_PYTHON_MODULE != 1
      }
    } // python
  } // LaDa
#endif

#ifdef LADA_INLINE
# undef LADA_INLINE
#endif
#ifdef LADA_END
# undef LADA_END
#endif
#if LADA_PYTHON_MODULE == 100
#  undef LADA_PYTHON_MODULE
#endif
// get ready for second inclusion
#ifdef LADA_PYTHON_MODULE 
# if LADA_PYTHON_MODULE == 0
#   undef LADA_PYTHON_MODULE 
#   define LADA_PYTHON_MODULE 1
# elif LADA_PYTHON_MODULE == 1
#   undef LADA_PYTHON_MODULE 
#   define LADA_PYTHON_MODULE 0
# endif
#endif 

#endif
