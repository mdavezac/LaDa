#ifndef PYLADA_CRYSTALMODULE_H
#define PYLADA_CRYSTALMODULE_H
#ifndef __cplusplus
# error LaDa requires a cpp compiler
#endif


#ifndef LADA_CRYSTAL_MODULE
#  define LADA_CRYSTAL_MODULE 100
#endif 

#if LADA_CRYSTAL_MODULE != 1
# include "LaDaConfig.h"
# include <Python.h>

# include <vector>
# include <cmath>

# include <boost/type_traits/is_same.hpp>
# include <boost/static_assert.hpp>
# include <boost/preprocessor/arithmetic/inc.hpp>
# include <python/ppslot.hpp>
# define BOOST_PP_VALUE 0
# include LADA_ASSIGN_SLOT(crystal)

# include <Eigen/LU> 

# include <root_exceptions.h>
# include <python/python.h>
# include <math/math.h>

# include <errors/exceptions.h>

# if LADA_CRYSTAL_MODULE == 100
    namespace LaDa
    {

      namespace crystal
      {
        /* This section is used in modules that use lada.crystal's API */
        static void **api_capsule = NULL; 
        
        namespace 
        {
          // Return -1 on error, 0 on success.
          // PyCapsule_Import will set an exception if there's an error.
          inline bool import(void)
          {
//           PyObject* module = PyImport_ImportModule("lada.crystal.cppwrappers");
//           if(not module) return false;
//           Py_DECREF(module);
            api_capsule = (void **)PyCapsule_Import("lada.crystal.cppwrappers._C_API", 0);
            return api_capsule != NULL;
          }
        }
      }
    }
# endif
#else
# define BOOST_PP_VALUE 0
# include LADA_ASSIGN_SLOT(crystal)
#endif

#if LADA_CRYSTAL_MODULE != 1
#  ifdef LADA_INLINE
#    error LADA_INLINE already defined
#  endif
#  if LADA_CRYSTAL_MODULE == 100
#    define LADA_INLINE inline
#  elif LADA_CRYSTAL_MODULE == 0
#    define LADA_INLINE
#  endif
#  ifdef LADA_END
#    error LADA_END already defined
#  elif LADA_CRYSTAL_MODULE == 0
#    define LADA_END(X) ;
#  elif LADA_CRYSTAL_MODULE == 100
#    define LADA_END(X) X
#  endif
#endif

#if LADA_CRYSTAL_MODULE != 1
  namespace LaDa
  {
    namespace crystal 
    {
      namespace 
      {
#endif

#       include "atom/atom.h"
#       include "structure/structure.h"
#       include "hart-forcade/hart-forcade.h"
#       include "utilities.h"
#       include "supercell.h"
#       include "algorithms.h"
        
#if LADA_CRYSTAL_MODULE != 1
      }
    }
  }
#endif


// get ready for second inclusion
#ifdef LADA_CRYSTAL_MODULE 
# if LADA_CRYSTAL_MODULE == 0
#   undef LADA_CRYSTAL_MODULE 
#   define LADA_CRYSTAL_MODULE 1
# elif LADA_CRYSTAL_MODULE == 1
#   undef LADA_CRYSTAL_MODULE 
#   define LADA_CRYSTAL_MODULE 0
# endif
#endif 
#ifdef LADA_INLINE
# undef LADA_INLINE
#endif
#ifdef LADA_END
# undef LADA_END
#endif

#if LADA_CRYSTAL_MODULE == 100
# undef LADA_CRYSTAL_MODULE
#endif

#endif 
