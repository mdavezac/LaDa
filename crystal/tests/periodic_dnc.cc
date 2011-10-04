#include "LaDaConfig.h"
#include <iostream>
#include <string>
#include <vector>
#include <set>

#include <cstdlib>
#include <time.h>

#include <boost/foreach.hpp>

#include <opt/debug.h>
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
                                LaDa::math::rVector3d const &_pos,
                                LaDa::math::iVector3d const &_n )
{
  using namespace LaDa::math;
  rVector3d const frac(_invcell * _pos);
  iVector3d const __ifrac(floor_int(frac));
  iVector3d const ifrac
    (
      __ifrac(0) + (__ifrac(0) < 0 ? _n(0): (__ifrac(0) >= _n(0)? -_n(0): 0)), 
      __ifrac(1) + (__ifrac(1) < 0 ? _n(1): (__ifrac(1) >= _n(1)? -_n(1): 0)), 
      __ifrac(2) + (__ifrac(2) < 0 ? _n(2): (__ifrac(2) >= _n(2)? -_n(2): 0))
    );
  iVector3d const neg(__ifrac(0) % _n(0), __ifrac(1) % _n(1), __ifrac(2) % _n(2));
  iVector3d const other
      ( 
        neg(0) < 0 ? neg(0) + _n(0): neg(0),
        neg(1) < 0 ? neg(1) + _n(1): neg(1),
        neg(2) < 0 ? neg(2) + _n(2): neg(2)
      );
// LADA_DOASSERT(other == ifrac, ~other << " | " << ~ifrac << " | " << ~_n << " | " << ~__ifrac << "\n");
  return other;
}

void check(LaDa::crystal::TemplateStructure< LADA_TYPE > const &_structure)
{
  using namespace LaDa;
  using namespace LaDa::math;
  using namespace LaDa::crystal;
  iVector3d const params = guess_mesh(_structure, 30);
  rMatrix3d invcell(gruber(_structure.cell()));
  for(size_t i(0); i < 3; ++i) invcell.col(i) /= types::t_real(params(i));
  invcell = invcell.inverse().eval();

  DnCBoxes dnc;
  dnc.init(_structure, params, 0.125);
  DnCBoxes::const_iterator i_box = dnc.begin();
  DnCBoxes::const_iterator i_box_end = dnc.end();
  std::cout << ~params << "\n";
  std::cout << _structure.cell() << "\n";
  for(size_t i(0); i_box != i_box_end; ++i_box, ++i)
  {
    std::cout << "i: " << i << "\n";
    DnCBoxes::value_type::const_iterator i_point = i_box->begin();
    DnCBoxes::value_type::const_iterator const i_point_end = i_box->end();
    for(; i_point != i_point_end; ++i_point)
      if(i_point->in_small_box) break;
    if(i_point == i_point_end)
    {
      std::cerr << "No points in box.\n";
      throw 0;
    }
    iVector3d const index = indices_(invcell, i_point->translation + _structure[i_point->index]->pos, params);
    for(i_point = i_box->begin(); i_point != i_point_end; ++i_point)
    {
      iVector3d const pi( indices_(invcell, i_point->translation + _structure[i_point->index]->pos, params) );
      if(i_point->in_small_box)
        { LADA_DOASSERT(math::eq(pi, index ), "Not in same box.\n"); }
      else
      {
        bool found = false;
        for(size_t i(0); i < 3; ++i)
          if(std::abs(pi(i) - index(i)) == 1) found = true;
          else if( params(i) > 1 and std::abs(pi(i) - index(i)) == params(i) - 1) found = true;
          else if( params(i) == 1 ) found = true;
        LADA_DOASSERT(found, "point in large box also in small box.\n")
      }
    }
  }
}

int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  using namespace LaDa::math;
  TemplateStructure< LADA_TYPE > lattice = b5(0.36); 

  check(lattice);


  types::t_int const seed = 1317585876; //time(NULL); // error 1317585909; // long 1317585876
  std::cout << seed << "\n";
  srand(seed);
  math::rMatrix3d cell;
  for(size_t k(0); k < 10; ++k)
  {
    for(size_t i(0); i < 3; ++i)
      for(size_t j(0); j < 3; ++j)
        cell(i, j) = rand() % 20 - 10;
    TemplateStructure< LADA_TYPE > structure = supercell(lattice, lattice.cell() * cell);
    std::cout << "natoms: " << structure.size() << "\n";
    check(structure);
  }

  return 0;
}
