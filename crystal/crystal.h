#ifndef PYLADA_CRYSTALMODULE_H
#define PYLADA_CRYSTALMODULE_H
#ifndef __cplusplus
# error LaDa requires a cpp compiler
#endif

#include <Python.h>

#ifndef LADA_CRYSTAL_MODULE
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

# include <boost/preprocessor/arithmetic/inc.hpp>
# include <boost/preprocessor/slot/slot.hpp>
# define BOOST_PP_VALUE 0
# include BOOST_PP_ASSIGN_SLOT(1)
#else
# include <boost/preprocessor/list/append.hpp>
# ifdef LADA_PYLIST
#   error LADA_PYLIST already defined
# endif
# ifdef LADA_PYLISTSAVE
#   error LADA_PYLISTSAVE already defined
# endif
# ifdef LADA_VALUE
#   error LADA_VALUE already defined
# endif
# define LADA_PYLIST (BOOST_PP_NIL)
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

#endif 
