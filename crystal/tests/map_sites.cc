#include "LaDaConfig.h"
#include <iostream>
#include <string>
#include <vector>
#include <set>

#include <opt/debug.h>
#include "../map_sites.h"
#include "../supercell.h"

namespace LaDa
{
  namespace crystal
  {
    TemplateStructure< LADA_TYPE > b5(types::t_real u)
    {
      TemplateStructure< LADA_TYPE > lattice;
      types::t_real const x(u), y(0.25 - u);
      lattice.set_cell(0, 0.5, 0.5)
                      (0.5, 0, 0.5)
                      (0.5, 0.5, 0);
      lattice.add_atom(5.000000e-01, 5.000000e-01, 5.000000e-01, "A")
                      (5.000000e-01, 2.500000e-01, 2.500000e-01, "A")
                      (2.500000e-01, 5.000000e-01, 2.500000e-01, "A")
                      (2.500000e-01, 2.500000e-01, 5.000000e-01, "A")
                      (8.750000e-01, 8.750000e-01, 8.750000e-01, "B")
                      (1.250000e-01, 1.250000e-01, 1.250000e-01, "B")
                      (     x,     x,     x, "X")
                      (     x,     y,     y, "X")
                      (     y,     x,     y, "X")
                      (     y,     y,     x, "X")
                      (    -x,    -x,    -x, "X")
                      (    -x,    -y,    -y, "X")
                      (    -y,    -x,    -y, "X")
                      (    -y,    -y,    -x, "X");
      return lattice;
    }
  }
}

void set(std::string &_in, std::string const &_type)
  { _in = _type; }
void set(std::vector<std::string> &_in, std::string const &_type)
  { _in[0] = _type; }
void set(std::set<std::string> &_in, std::string const &_type)
  { _in.clear(); _in.insert(_type); }

int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  using namespace LaDa::math;
  TemplateStructure< LADA_TYPE > lattice = b5(0.36); 
  math::rMatrix3d cell;
  cell << 5, 2, -5, 0, 2, 3, -3, 0, 1;
  lattice.scale() = 2e0;
  TemplateStructure< LADA_TYPE > super = supercell(lattice, lattice.cell() * cell);
  TemplateStructure< LADA_TYPE > checkme = super.copy();
  checkme.scale() *= 0.5;
  checkme.cell() *= 2e0;
  for(size_t i(0); i < checkme.size(); ++i) checkme[i]->pos *= 2e0;;

  LADA_DOASSERT(map_sites(lattice, checkme) == true, "Should have mapped all sites.\n"); 
  for(size_t i(0); i < checkme.size(); ++i)
    LADA_DOASSERT(checkme[i]->site == super[i]->site, "Site indices do not correspond.");

  set(checkme.front()->type, "Ge");
  LADA_DOASSERT(map_sites(lattice, checkme) == false, "Should not have mapped all sites."); 
  for(size_t i(1); i < checkme.size()-1; ++i)
    LADA_DOASSERT(checkme[i]->site == super[i]->site, "Site indices do not correspond.");
  LADA_DOASSERT(checkme[0]->site == -1, "Site should not be mapped.");
 
  checkme.add_atom(0.5, 0.5, 0.125, "Si"); 
  LADA_DOASSERT(map_sites(lattice, checkme, false) == false, "Should not have mapped all sites."); 
  for(size_t i(0); i < checkme.size()-1; ++i)
    LADA_DOASSERT(checkme[i]->site == super[i]->site, "Site indices do not correspond.");
  LADA_DOASSERT(checkme[checkme.size()-1]->site == -1, "Site should not be mapped.");

  return 0;
}
