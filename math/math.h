#ifndef PYLADA_MATH_MODULE_H
#define PYLADA_MATH_MODULE_H
#ifndef __cplusplus
# error LaDa requires a cpp compiler
#endif


#ifndef LADA_MATH_MODULE
#  define LADA_MATH_MODULE 100
#endif 

#if LADA_MATH_MODULE != 1
# include "LaDaConfig.h"
# include <Python.h>
# ifndef LADA_PYTHONTWOSIX
#   if PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7
#     define LADA_PYTHONTWOSIX
#   endif
# endif

# include <boost/utility/enable_if.hpp>
# include <boost/type_traits/is_integral.hpp>
# include <boost/type_traits/is_arithmetic.hpp>
# include <boost/type_traits/is_floating_point.hpp>
# include <boost/numeric/conversion/converter.hpp>

# include <boost/preprocessor/arithmetic/inc.hpp>
# include <python/ppslot.hpp>
# define BOOST_PP_VALUE 0
# include LADA_ASSIGN_SLOT(math)

# include <python/types.h>
# include <Eigen/Core>
# include <Eigen/LU>
# include <Eigen/Geometry>

# include <errors/exceptions.h>

# if LADA_MATH_MODULE == 100

    namespace LaDa
    {

      namespace math
      {
        /* This section is used in modules that use lada.math's API */
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
            PyObject *module = PyImport_ImportModule("lada.math");
            if(not module) return false;
#           ifdef LADA_PYTHONTWOSIX
              PyObject* c_api_object = PyObject_GetAttrString(module, "_C_API");
	      if (c_api_object == NULL) { Py_DECREF(module); return false; }
              if (PyCObject_Check(c_api_object))
                api_capsule = (void **)PyCObject_AsVoidPtr(c_api_object);
              Py_DECREF(c_api_object);
#           else
              api_capsule = (void **)PyCapsule_Import("lada.math._C_API", 0);
#           endif
            Py_DECREF(module);
            return api_capsule != NULL;
          }
        }
      }
    }
# endif
#else
# define BOOST_PP_VALUE 0
# include LADA_ASSIGN_SLOT(math)
#endif

#if LADA_MATH_MODULE != 1
#  ifdef LADA_INLINE
#    error LADA_INLINE already defined
#  endif
#  if LADA_MATH_MODULE == 100
#    define LADA_INLINE inline
#  elif LADA_MATH_MODULE == 0
#    define LADA_INLINE
#  endif
#  ifdef LADA_END
#    error LADA_END already defined
#  elif LADA_MATH_MODULE == 0
#    define LADA_END(X) ;
#  elif LADA_MATH_MODULE == 100
#    define LADA_END(X) { X }
#  endif
#endif

// Some stuff needs to be added to eigen, so do outside of namespace.
// Also, it is actually needed by python. Not very cool. 
# include "eigen.h"
// different namespace
# include "exceptions.h"

#if LADA_MATH_MODULE != 1
  namespace LaDa
  {
    namespace math 
    {
      namespace 
      {
#endif

#       include "fuzzy.h"
#       include "misc.h"
#       include "symmetry_operator.h"
#       include "algorithms.h"
        
#if LADA_MATH_MODULE != 1
      }
    }
  }
#endif


// get ready for second inclusion
#ifdef LADA_MATH_MODULE 
# if LADA_MATH_MODULE == 0
#   undef LADA_MATH_MODULE 
#   define LADA_MATH_MODULE 1
# elif LADA_MATH_MODULE == 1
#   undef LADA_MATH_MODULE 
#   define LADA_MATH_MODULE 0
# endif
#endif 
#ifdef LADA_INLINE
# undef LADA_INLINE
#endif
#ifdef LADA_END
# undef LADA_END
#endif

#if LADA_MATH_MODULE == 100
# undef LADA_MATH_MODULE
#endif

#endif 

