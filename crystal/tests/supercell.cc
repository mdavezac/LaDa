#include "LaDaConfig.h"
#include<iostream>
#include<string>
#include<vector>

#include <opt/debug.h>
#include "../supercell.h"
#include "../compare_sites.h"

#if LADA_INCREMENT == 0
#  define LADA_INIT 
#  define LADA_DOASSERT_TYPE0 LADA_DOASSERT(i_atom->type() == "Si", "Wrong first type.\n"); 
#  define LADA_DOASSERT_TYPE1 LADA_DOASSERT(i_atom->type() == "Ge", "Wrong first type.\n"); 
#elif LADA_INCREMENT == 1
#  define LADA_INIT , "Si"
#  define LADA_DOASSERT_TYPE0 \
     LADA_DOASSERT(i_atom->type().size() == 1, "Wrong number of types.\n");\
     LADA_DOASSERT(i_atom->type()[0] == "Si", "Wrong first type.\n");
#  define LADA_DOASSERT_TYPE1 \
     LADA_DOASSERT(i_atom->type().size() == 2, "Wrong number of types.\n");\
     LADA_DOASSERT(i_atom->type()[0] == "Ge", "Wrong first type.\n");\
     LADA_DOASSERT(i_atom->type()[1] == "Si", "Wrong first type.\n");
#elif LADA_INCREMENT == 2
#  define LADA_INIT , "Si"
#  define LADA_DOASSERT_TYPE0 \
    { LADA_TYPE cmp; cmp.insert("Si"); \
      LADA_DOASSERT(i_atom->type().size() == 1, "Wrong number of types.\n");\
      LADA_DOASSERT(i_atom->type() == cmp, "Wrong first type.\n"); }
#  define LADA_DOASSERT_TYPE1 \
    { LADA_TYPE cmp; cmp.insert("Si"); cmp.insert("Ge");                \
      LADA_DOASSERT(i_atom->type().size() == 2, "Wrong number of types.\n");\
      LADA_DOASSERT(i_atom->type() == cmp, "Wrong first type.\n"); }
#endif

using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  using namespace LaDa::math;
  Structure< LADA_TYPE > lattice;
  lattice.set_cell(0,0.5,0.5)
                  (0.5,0,0.5)
                  (0.5,0.5,0);
  lattice.add_atom(0,0,0, "Si") 
                  (0.25,0.25,0.25, "Ge" LADA_INIT);
  lattice[0]->freeze = AtomFreezeMixin::frozen::Y;

  rMatrix3d matrix;
  matrix << -1, 1, 1, 1, -1, 1, 1, 1, -1;

  Structure< LADA_TYPE > result;
  result = supercell(lattice, lattice.cell() * matrix);
  LADA_DOASSERT(is_identity(result.cell()), "Unexpected supercell.");
  LADA_DOASSERT(eq(result[0]->pos, rVector3d(0.00000, 0.00000, 0.00000)), "Incorrect position.\n")
  LADA_DOASSERT(eq(result[1]->pos, rVector3d(0.25000, 0.25000, 0.25000)), "Incorrect position.\n")
  LADA_DOASSERT(eq(result[2]->pos, rVector3d(0.50000, 0.00000, 0.50000)), "Incorrect position.\n")
  LADA_DOASSERT(eq(result[3]->pos, rVector3d(0.75000, 0.25000, 0.75000)), "Incorrect position.\n")
  LADA_DOASSERT(eq(result[4]->pos, rVector3d(0.50000, 0.50000, 0.00000)), "Incorrect position.\n")
  LADA_DOASSERT(eq(result[5]->pos, rVector3d(0.75000, 0.75000, 0.25000)), "Incorrect position.\n")
  LADA_DOASSERT(eq(result[6]->pos, rVector3d(0.00000, 0.50000, 0.50000)), "Incorrect position.\n")
  LADA_DOASSERT(eq(result[7]->pos, rVector3d(0.25000, 0.75000, 0.75000)), "Incorrect position.\n")
  Structure< LADA_TYPE > :: const_iterator i_atom = result.begin();
  Structure< LADA_TYPE > :: const_iterator const i_atom_end = result.begin();
  for(; i_atom != i_atom_end; ++i_atom)
  {
    if(are_periodic_images(lattice[0]->pos, i_atom->pos(), lattice.cell().inverse()))
    {
      LADA_DOASSERT(i_atom->freeze() == AtomFreezeMixin::frozen::Y, "Incorrect freeze.\n");
      LADA_DOASSERT(i_atom->site() == 0, "Incorrect site index.\n"); 
      LADA_DOASSERT_TYPE0;
    }
    else if(are_periodic_images(lattice[1]->pos, i_atom->pos(), lattice.cell().inverse()))
    {
      LADA_DOASSERT(i_atom->freeze() == AtomFreezeMixin::frozen::NONE, "Incorrect freeze.\n");
      LADA_DOASSERT(i_atom->site() == 1, "Incorrect site index.\n"); 
      LADA_DOASSERT_TYPE1;
    }
    else { LADA_DOASSERT(true, "Incorrect site index.\n"); }
  }
  return 0;
}
