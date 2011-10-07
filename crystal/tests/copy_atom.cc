#include "LaDaConfig.h"
#include <iostream>
#include <string>
#include <vector>
#include <set>

#include <opt/debug.h>
#include "../atom.h"

#if LADA_INCREMENT == 0
#  define LADA_INIT_TYPE atom.type() = "hello"; atom->site = 2; 
#  define LADA_DOASSERT_TYPE LADA_DOASSERT(c.type() == "hello", "Wrong first type.\n");
#elif LADA_INCREMENT == 1
#  define LADA_INIT_TYPE atom.type().push_back("Au"); atom.type().push_back("Pd"); atom->site = 2;
#  define LADA_DOASSERT_TYPE \
     LADA_DOASSERT(c.type().size() == 2, "Wrong number of types.\n");\
     LADA_DOASSERT(c.type()[0] == "Au", "Wrong first type.\n");\
     LADA_DOASSERT(c.type()[1] == "Pd", "Wrong first type.\n");
#elif LADA_INCREMENT == 2
#  define LADA_INIT_TYPE atom.type().insert("Au"); atom.type().insert("Pd"); atom->site = 2;
#  define LADA_DOASSERT_TYPE \
    { LADA_TYPE cmp; cmp.insert("Au"); cmp.insert("Pd");                \
      LADA_DOASSERT(c.type().size() == 2, "Wrong number of types.\n");\
      LADA_DOASSERT(c.type() == cmp, "Wrong first type.\n"); }
#endif
using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  Atom< LADA_TYPE > atom;
  atom.pos() = math::rVector3d(0.5,-0.5,0.5);
  atom.freeze() = AtomFreezeMixin::frozen::X;
  LADA_INIT_TYPE;

  Atom< LADA_TYPE > c = atom.copy();
  LADA_DOASSERT_TYPE;
  LADA_DOASSERT(c->site == 2, "Wrong site index.\n");
  LADA_DOASSERT( (c->pos-atom->pos).squaredNorm() < 1e-12, "Wrong positions.\n");

  return 0;
}
