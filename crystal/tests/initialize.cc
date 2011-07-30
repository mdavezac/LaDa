

#include "LaDaConfig.h"
#include<iostream>
#include<string>
#include<vector>
#include "../structure.h"
#include "../is_container.h"
#include <opt/debug.h>

using namespace std;
int main()
{
  using namespace LaDa::crystal;
# if LADA_TEST_INCTYPE == 0
    LADA_TEST_STRUCTURE<std::string> structure;
# elif LADA_TEST_INCTYPE == 1
    LADA_TEST_STRUCTURE< std::vector<std::string> > structure;
# elif LADA_TEST_INCTYPE == 2
    LADA_TEST_STRUCTURE<int> structure;
# endif
  structure.set_cell(-0.5,0.5,0.5)
                    (0.5,-0.5,0.5)
                    (0.5,0.5,-0.5);
# if LADA_TEST_INCTYPE == 0
    structure.add_atom(0,0,0, "Au");
    structure.add_atom(0.25,0.25,0.25, "Pd");
# elif LADA_TEST_INCTYPE == 1
    structure.add_atom(0,0,0, "Au");
    structure.add_atom(0.25,0.25,0.25, "Au", "Pd");
# elif LADA_TEST_INCTYPE == 2
    structure.add_atom(0,0,0, 1);
    structure.add_atom(0.25,0.25,0.25, 15);
# endif
  for(size_t i(0); i < 3; ++i)
  {
    for(size_t j(0); j < 3; ++j)
      if(i == j) { LADA_DOASSERT(structure(i,j) == -0.5, "structure cell data incorrect."); }
      else { LADA_DOASSERT(structure(i,j) == 0.5, "structure cell data incorrect."); }
    LADA_DOASSERT(structure[0].pos[i] == 0, "structure atom data incorrect.");
    LADA_DOASSERT(structure[1].pos[i] == 0.25, "structure atom data incorrect.");
  }
# if LADA_TEST_INCTYPE == 0
    LADA_DOASSERT(structure[0].type == "Au", "structure atom data incorrect.");
    LADA_DOASSERT(structure[1].type == "Pd", "structure atom data incorrect.");
# elif LADA_TEST_INCTYPE == 1
    LADA_DOASSERT(structure[0].type.size() == 1, "structure atom data incorrect.");
    LADA_DOASSERT(structure[0].type[0] == "Au", "structure atom data incorrect.");
    LADA_DOASSERT(structure[1].type.size() == 2, "structure atom data incorrect.");
    LADA_DOASSERT(structure[1].type[0] == "Au", "structure atom data incorrect.");
    LADA_DOASSERT(structure[1].type[1] == "Pd", "structure atom data incorrect.");
# elif LADA_TEST_INCTYPE == 2
    LADA_DOASSERT(structure[0].type == 1, "structure atom data incorrect.");
    LADA_DOASSERT(structure[1].type == 15, "structure atom data incorrect.");
# endif
  return 0;
}
