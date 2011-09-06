#include "LaDaConfig.h"
#include <iostream>
#include <string>
#include <vector>
#include <set>

#include <boost/foreach.hpp>

#include "../periodic_dnc.h"
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

void set(std::vector<std::string> &_in, std::string const &_type)
  { _in[0] = _type; }
void set(std::set<std::string> &_in, std::string const &_type)
  { _in.clear(); _in.insert(_type); }

LaDa::math::iVector3d indices_( LaDa::math::rMatrix3d const &_invcell,
                                LaDa::math::rVector3d const &_pos )
{
  using namespace LaDa::math;
  rVector3d const frac(_invcell * _pos);
  return iVector3d(std::floor(frac(0)+1e-12), std::floor(frac(1)+1e-12), std::floor(frac(1)+1e-12));
}

int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  using namespace LaDa::math;
  TemplateStructure< LADA_TYPE > lattice = b5(0.36); 
  math::rMatrix3d cell;
  cell << 5, 0, 0, 0, 1, 0, 0, 0, 1;
  TemplateStructure< LADA_TYPE > structure = supercell(lattice, lattice.cell() * cell);
  cell(0,0) = 1;
  math::rMatrix3d invcell(cell.inverse());
  DnCBoxes dnc;
  dnc.init(structure, iVector3d(5, 1, 1), 0.125);
  DnCBoxes::const_iterator i_box = dnc.begin();
  DnCBoxes::const_iterator i_box_end = dnc.end();
  for(size_t i(0); i_box != i_box_end; ++i_box, ++i)
  {
    DnCBoxes::value_type::const_iterator i_point = i_box->begin();
    DnCBoxes::value_type::const_iterator const i_point_end = i_box->end();
    for(; i_point != i_point_end; ++i_point)
      if(i_point->in_small_box) break;
    LADA_DOASSERT(i_point != i_point_end, "No points in box.\n")
    iVector3d indices = indices_(invcell, i_point->translation + structure[i_point->index].pos);
    std::cout << ~indices << "\n";
    continue;
    for(i_point = i_box->begin(); i_point != i_point_end; ++i_point)
      if(i_point->in_small_box)
      {
        LADA_DOASSERT(math::eq( indices_(invcell, i_point->translation + structure[i_point->index].pos),
                                indices ), "Not in same box.\n");
      }
      else
      {
        iVector3d const other = indices_(invcell, i_point->translation + structure[i_point->index].pos);
        std::cout << "   - " << ~other << "\n";
        LADA_DOASSERT(    std::abs(other(0) - indices(0)) == 5 
                       or std::abs(other(0) - indices(0)) == 1, "Found in unexpected box.\n");
        LADA_DOASSERT(other(1) == indices(1), "Found in unexpected box.\n");
        LADA_DOASSERT(other(2) == indices(1), "Found in unexpected box.\n");
      }
  }

  return 0;
}
