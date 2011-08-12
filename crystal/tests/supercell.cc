#include "LaDaConfig.h"
#include<iostream>
#include<string>
#include<vector>

#include <boost/foreach.hpp>
#include "../supercell.h"

#define TYPE std::vector<std::string>

using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  using namespace LaDa::math;
  TemplateStructure< TYPE > lattice;
  lattice.set_cell(0,0.5,0.5)
                  (0.5,0,0.5)
                  (0.5,0.5,0);
  lattice.add_atom(0,0,0, "Ge") 
         .add_atom(0.25,0.25,0.25, "Si","Ge");
  lattice[0].freeze = Atom<TYPE>::frozen::X;

  rMatrix3d matrix;
  matrix << -1, 1, 1, 1, -1, 1, 1, 1, -1;

  TemplateStructure< std::vector<std::string> > result;
  result = supercell(lattice, lattice.cell() * matrix);
  LADA_ASSERT(is_identity(result.cell()), "Unexpected supercell.");
  LADA_ASSERT(eq(result[0].pos, rVector3d(0.00000, 0.00000, 0.00000)), "Incorrect position.\n")
  LADA_ASSERT(eq(result[1].pos, rVector3d(0.25000, 0.25000, 0.25000)), "Incorrect position.\n")
  LADA_ASSERT(eq(result[2].pos, rVector3d(0.50000, 0.00000, 0.50000)), "Incorrect position.\n")
  LADA_ASSERT(eq(result[3].pos, rVector3d(0.75000, 0.25000, 0.75000)), "Incorrect position.\n")
  LADA_ASSERT(eq(result[4].pos, rVector3d(0.50000, 0.50000, 0.00000)), "Incorrect position.\n")
  LADA_ASSERT(eq(result[5].pos, rVector3d(0.75000, 0.75000, 0.25000)), "Incorrect position.\n")
  LADA_ASSERT(eq(result[6].pos, rVector3d(0.00000, 0.50000, 0.50000)), "Incorrect position.\n")
  LADA_ASSERT(eq(result[7].pos, rVector3d(0.25000, 0.75000, 0.75000)), "Incorrect position.\n")
  foreach(Atom<TYPE> const &atom, result)
  {
    if(are_periodic_images(lattice[0].pos, atom.pos, lattice.cell().inverse()))
    {
      LADA_ASSERT(atom.freeze == Atom<TYPE>::frozen::X, "Incorrect freeze.\n");
      LADA_ASSERT(atom.site == 0, "Incorrect site index.\n"); 
      LADA_ASSERT(atom.type.size() == 1,  "Incorrect number of species.\n");
      LADA_ASSERT(atom.type[0] == "Ge",  "Incorrect specie.\n");
    }
    else if(are_periodic_images(lattice[1].pos, atom.pos, lattice.cell().inverse()))
    {
      LADA_ASSERT(atom.freeze == Atom<TYPE>::frozen::NONE, "Incorrect freeze.\n");
      LADA_ASSERT(atom.site == 1, "Incorrect site index.\n"); 
      LADA_ASSERT(atom.type.size() == 2,  "Incorrect number of species.\n");
      LADA_ASSERT(atom.type[0] == "Si",  "Incorrect specie.\n");
      LADA_ASSERT(atom.type[1] == "Ge",  "Incorrect specie.\n");
    }
    else { LADA_ASSERT(true, "Incorrect site index.\n"); }
  }
  return 0;
}
