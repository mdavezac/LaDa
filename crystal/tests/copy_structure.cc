#include "LaDaConfig.h"
#include <iostream>
#include <string>
#include <vector>
#include <set>

#include <opt/debug.h>
#include <math/misc.h>
#include "../structure.h"

#if LADA_INCREMENT == 0
#  define LADA_INIT 
#  define LADA_DOASSERT_TYPE0 LADA_DOASSERT(c[0]->type == "Si", "Wrong first type.\n"); 
#  define LADA_DOASSERT_TYPE1 LADA_DOASSERT(c[1]->type == "Ge", "Wrong first type.\n"); 
#elif LADA_INCREMENT == 1
#  define LADA_INIT , "Si"
#  define LADA_DOASSERT_TYPE0 \
     LADA_DOASSERT(c[0]->type.size() == 1, "Wrong number of types.\n");\
     LADA_DOASSERT(c[0]->type[0] == "Si", "Wrong first type.\n");
#  define LADA_DOASSERT_TYPE1 \
     LADA_DOASSERT(c[1]->type.size() == 2, "Wrong number of types.\n");\
     LADA_DOASSERT(c[1]->type[0] == "Ge", "Wrong first type.\n");\
     LADA_DOASSERT(c[1]->type[1] == "Si", "Wrong first type.\n");
#elif LADA_INCREMENT == 2
#  define LADA_INIT , "Si"
#  define LADA_DOASSERT_TYPE0 \
    { LADA_TYPE cmp; cmp.insert("Si"); \
      LADA_DOASSERT(c[0]->type.size() == 1, "Wrong number of types.\n");\
      LADA_DOASSERT(c[0]->type == cmp, "Wrong first type.\n"); }
#  define LADA_DOASSERT_TYPE1 \
    { LADA_TYPE cmp; cmp.insert("Si"); cmp.insert("Ge");                \
      LADA_DOASSERT(c[1]->type.size() == 2, "Wrong number of types.\n");\
      LADA_DOASSERT(c[1]->type == cmp, "Wrong first type.\n"); }
#endif
using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  using namespace LaDa::math;
  TemplateStructure< LADA_TYPE > lattice;
  lattice.set_cell(0,0.5,0.5)
                  (0.5,0,0.5)
                  (0.5,0.5,0);
  lattice.add_atom(0,0,0, "Si")
                  (0.25,0.25,0.25, "Ge" LADA_INIT);
  lattice[0]->freeze = AtomFreezeMixin::frozen::X;
  lattice[1]->freeze = AtomFreezeMixin::frozen::Y;
  lattice[0]->site = 0;
  lattice[1]->site = 1;
  lattice->energy = -1;
  lattice->weight = -1;
  lattice->scale = 0.5;
  lattice->freeze = frozenstr::XX;
  lattice->name = "hello";

  TemplateStructure< LADA_TYPE > c = lattice.copy();
  LADA_DOASSERT(is_null(c->cell - lattice->cell), "Different cells.\n");
  LADA_DOASSERT(is_null(c->energy - lattice->energy), "Different energies.\n");
  LADA_DOASSERT(is_null(c->weight - lattice->weight), "Different weights.\n");
  LADA_DOASSERT(is_null(c->scale - lattice->scale), "Different scale.\n");
  LADA_DOASSERT(c->name == lattice->name, "Different names.\n");
  LADA_DOASSERT(c->freeze == lattice->freeze, "Different frozen dofs.\n");
  LADA_DOASSERT(c.size() == lattice.size(), "Different sizes.\n");
  LADA_DOASSERT(is_null(c[0]->pos - lattice[0]->pos), "Different atomic position 0.\n");
  LADA_DOASSERT(is_null(c[0]->pos - lattice[0]->pos), "Different atomic position 0.\n");
  LADA_DOASSERT(c[0]->freeze == lattice[0]->freeze, "Different frozen dofs at site 0.\n");
  LADA_DOASSERT(c[1]->freeze == lattice[1]->freeze, "Different frozen dofs at site 1.\n");
  LADA_DOASSERT(c[0]->site == lattice[0]->site, "Different site index at site 0.\n");
  LADA_DOASSERT(c[1]->site == lattice[1]->site, "Different site index at site 0.\n");
  LADA_DOASSERT_TYPE0;
  LADA_DOASSERT_TYPE1;

  return 0;
}
