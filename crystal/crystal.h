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
      
      // Return -1 on error, 0 on success.
      // PyCapsule_Import will set an exception if there's an error.
      static bool import(void)
      {
        api_capsule = (void **)PyCapsule_Import("lada.crystal.cppwrappers._C_API", 0);
        return api_capsule != NULL;
      }
    }
  }
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
