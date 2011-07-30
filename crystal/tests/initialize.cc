

#include "LaDaConfig.h"
#include<iostream>
#include<string>
#include<vector>
#include "../structure.h"
#include "../is_container.h"
#include <opt/debug.h>

void check_structuredata()
{
  using namespace LaDa::crystal;
  StructureData<std::string> structure;
  structure.set_cell(-0.5,0.5,0.5)
               (0.5,-0.5,0.5)
               (0.5,0.5,-0.5);
  structure.add_atom(0,0,0, "Au");
  structure.add_atom(0.25,0.25,0.25, "Pd");
  for(size_t i(0); i < 3; ++i)
  {
    for(size_t j(0); j < 3; ++j)
      if(i == j) { LADA_DOASSERT(structure.cell(i,j) == -0.5, "structure cell data incorrect."); }
      else { LADA_DOASSERT(structure.cell(i,j) == 0.5, "structure cell data incorrect."); }
    LADA_DOASSERT(structure[0].pos[i] == 0, "structure atom data incorrect.");
    LADA_DOASSERT(structure[1].pos[i] == 0.25, "structure atom data incorrect.");
  }
  LADA_DOASSERT(structure[0].type == "Au", "structure atom data incorrect.");
  LADA_DOASSERT(structure[1].type == "Pd", "structure atom data incorrect.");
}
void check_structure()
{
  using namespace LaDa::crystal;
  TemplateStructure<std::string> structure;
  structure.set_cell(-0.5,0.5,0.5)
               (0.5,-0.5,0.5)
               (0.5,0.5,-0.5);
  structure.add_atom(0,0,0, "Au");
  structure.add_atom(0.25,0.25,0.25, "Pd");
  for(size_t i(0); i < 3; ++i)
  {
    for(size_t j(0); j < 3; ++j)
      if(i == j) { LADA_DOASSERT(structure(i,j) == -0.5, "structure cell data incorrect."); }
      else { LADA_DOASSERT(structure(i,j) == 0.5, "structure cell data incorrect."); }
    LADA_DOASSERT(structure[0].pos[i] == 0, "structure atom data incorrect.");
    LADA_DOASSERT(structure[1].pos[i] == 0.25, "structure atom data incorrect.");
  }
  LADA_DOASSERT(structure[0].type == "Au", "structure atom data incorrect.");
  LADA_DOASSERT(structure[1].type == "Pd", "structure atom data incorrect.");
}

void check_vstructuredata()
{
  using namespace LaDa::crystal;
  StructureData< std::vector<std::string> > structure;
  structure.set_cell(-0.5,0.5,0.5)
               (0.5,-0.5,0.5)
               (0.5,0.5,-0.5);
  structure.add_atom(0,0,0, "Au");
  structure.add_atom(0.25,0.25,0.25, "Au", "Pd");
  for(size_t i(0); i < 3; ++i)
  {
    for(size_t j(0); j < 3; ++j)
      if(i == j) { LADA_DOASSERT(structure.cell(i,j) == -0.5, "structure cell data incorrect."); }
      else { LADA_DOASSERT(structure.cell(i,j) == 0.5, "structure cell data incorrect."); }
    LADA_DOASSERT(structure[0].pos[i] == 0, "structure atom data incorrect.");
    LADA_DOASSERT(structure[1].pos[i] == 0.25, "structure atom data incorrect.");
  }
  LADA_DOASSERT(structure[0].type.size() == 1, "structure atom data incorrect.");
  LADA_DOASSERT(structure[0].type[0] == "Au", "structure atom data incorrect.");
  LADA_DOASSERT(structure[1].type.size() == 2, "structure atom data incorrect.");
  LADA_DOASSERT(structure[1].type[0] == "Au", "structure atom data incorrect.");
  LADA_DOASSERT(structure[1].type[1] == "Pd", "structure atom data incorrect.");
}

void check_vstructure()
{
  using namespace LaDa::crystal;
  TemplateStructure< std::vector<std::string> > structure;
  structure.set_cell(-0.5,0.5,0.5)
               (0.5,-0.5,0.5)
               (0.5,0.5,-0.5);
  structure.add_atom(0,0,0, "Au");
  structure.add_atom(0.25,0.25,0.25, "Au", "Pd");
  for(size_t i(0); i < 3; ++i)
  {
    for(size_t j(0); j < 3; ++j)
      if(i == j) { LADA_DOASSERT(structure(i,j) == -0.5, "structure cell data incorrect."); }
      else { LADA_DOASSERT(structure(i,j) == 0.5, "structure cell data incorrect."); }
    LADA_DOASSERT(structure[0].pos[i] == 0, "structure atom data incorrect.");
    LADA_DOASSERT(structure[1].pos[i] == 0.25, "structure atom data incorrect.");
  }
  LADA_DOASSERT(structure[0].type.size() == 1, "structure atom data incorrect.");
  LADA_DOASSERT(structure[0].type[0] == "Au", "structure atom data incorrect.");
  LADA_DOASSERT(structure[1].type.size() == 2, "structure atom data incorrect.");
  LADA_DOASSERT(structure[1].type[0] == "Au", "structure atom data incorrect.");
  LADA_DOASSERT(structure[1].type[1] == "Pd", "structure atom data incorrect.");
}
int main()
{
  using namespace LaDa::crystal;
  check_structuredata();
  check_vstructuredata();
  check_structure();
  check_vstructure();

  return 0;
};
