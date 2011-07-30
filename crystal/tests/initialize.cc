

#include "LaDaConfig.h"
#include<iostream>
#include<string>
#include<vector>
#include "../structure.h"
#include "../is_container.h"
#include <opt/debug.h>

#if LADA_STRUCTURE == 0
#  define T_STR StructureData
#else
#  define T_STR TemplateStructure
#endif
#if LADA_TYPE == 0
#  define T_TYPE std::string
#else
#  define T_TYPE std::vector<std::strin>
#endif

int main()
{
  using namespace LaDa::crystal;
  T_STR<T_TYPE> structure;
  structure.set_cell(-0.5,0.5,0.5)
               (0.5,-0.5,0.5)
               (0.5,0.5,-0.5);
  structure.add_atom(0,0,0, "Au");
# if LADA_TYPE == 0
    structure.add_atom(0.25,0.25,0.25, "Pd");
# else
    structure.add_atom(0.25,0.25,0.25, "Au", "Pd");
# endif
  for(size_t i(0); i < 3; ++i)
  {
    for(size_t j(0); j < 3; ++j)
      if(i == j) { LADA_DOASSERT(structure(i,j) == -0.5, "structure cell data incorrect."); }
      else { LADA_DOASSERT(structure(i,j) == 0.5, "structure cell data incorrect."); }
    LADA_DOASSERT(structure[0].pos[i] == 0, "structure atom data incorrect.");
    LADA_DOASSERT(structure[1].pos[i] == 0.25, "structure atom data incorrect.");
  }
# if LADA_TYPE == 0
    LADA_DOASSERT(structure[0].type == "Au", "structure atom data incorrect.");
    LADA_DOASSERT(structure[1].type == "Pd", "structure atom data incorrect.");
# else
    LADA_DOASSERT(structure[0].type.size() == 1, "structure atom data incorrect.");
    LADA_DOASSERT(structure[0].type[0] == "Au", "structure atom data incorrect.");
    LADA_DOASSERT(structure[1].type.size() == 2, "structure atom data incorrect.");
    LADA_DOASSERT(structure[1].type[0] == "Au", "structure atom data incorrect.");
    LADA_DOASSERT(structure[1].type[1] == "Pd", "structure atom data incorrect.");
# endif
  return 0;
}
