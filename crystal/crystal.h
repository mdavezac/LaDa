#ifndef PYLADA_CRYSTALMODULE_H
#define PYLADA_CRYSTALMODULE_H
#ifndef __cplusplus
# error LaDa requires a cpp compiler
#endif


#ifndef LADA_CRYSTAL_MODULE
#  define LADA_CRYSTAL_MODULE 100
#endif 

#if LADA_CRYSTAL_MODULE != 1
# include <boost/preprocessor/arithmetic/inc.hpp>
# include <boost/preprocessor/slot/slot.hpp>
# define BOOST_PP_VALUE 0
# include BOOST_PP_ASSIGN_SLOT(1)

# include <Python.h>
  namespace LaDa
  {
    namespace crystal
    {
      /* This section is used in modules that use lada.crystal's API */
      static void **api_capsule;
      
      namespace 
      {
        // Return -1 on error, 0 on success.
        // PyCapsule_Import will set an exception if there's an error.
        inline bool import(void)
        {
          api_capsule = (void **)PyCapsule_Import("lada.crystal.cppwrappers._C_API", 0);
          return api_capsule != NULL;
        }
      }
    }
  }
#endif
#if LADA_CRYSTAL_MODULE == 0
# include <boost/preprocessor/arithmetic/inc.hpp>
# include <boost/preprocessor/slot/slot.hpp>
# define BOOST_PP_VALUE 0
# include BOOST_PP_ASSIGN_SLOT(1)
#elif LADA_CRYSTAL_MODULE == 1
# define BOOST_PP_VALUE 0
# include BOOST_PP_ASSIGN_SLOT(1)
#endif

#if LADA_CRYSTAL_MODULE != 1
#  ifdef LADA_END
#    error LADA_END already defined
#  elif LADA_CRYSTAL_MODULE == 0
#    define LADA_END(X) ;
#  elif LADA_CRYSTAL_MODULE == 100
#    define LADA_END(X) X
#  endif
#endif

#include "atom/pybase.h"
#include "atom/atom.h"
#include "structure/pybase.h"
#include "structure/structure.h"
#include "hart-forcade/pybase.h"
#include "hart-forcade/hart-forcade.h"

#include "utilities.h"
#include "supercell.h"
#include "map_sites.h"
#include "equivalent_structures.h"
#include "primitive.h"
#include "space_group.h"
#include "neighbors.h"
#include "coordination_shells.h"
#include "confsplit.h"
#include "periodic_dnc.h"

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
#ifdef LADA_END
# undef LADA_END
#endif

#if LADA_CRYSTAL_MODULE == 100
# undef LADA_CRYSTAL_MODULE
#endif

#endif 
