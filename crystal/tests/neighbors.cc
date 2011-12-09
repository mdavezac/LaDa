#include "LaDaConfig.h"

#include<iostream>
#include<string>
#include<vector>
#include<set>

#include <boost/foreach.hpp>

#include "../neighbors.h"
#include "../supercell.h"

using namespace std;
void check( LaDa::crystal::Neighbors::const_iterator &i_first, 
            LaDa::crystal::Neighbors::const_iterator &i_end ) 
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  using namespace LaDa::math;
  for(size_t i(0); i < 4; ++i, ++i_first)
  {
    LADA_DOASSERT(i_first != i_end, "Premature death.\n");
    LADA_DOASSERT(eq(i_first->distance, sqrt(3.0)*0.25), "Wrong distance.\n");
    for(size_t i(0); i < 3; ++i)
      LADA_DOASSERT(eq(abs(i_first->pos(i)), 0.25), "Wrong position.\n");
  }
  for(size_t i(0); i < 12; ++i, ++i_first)
  {
    LADA_DOASSERT(i_first != i_end, "Premature death.\n");
    LADA_DOASSERT(eq(i_first->distance, sqrt(2.0)*0.5), "Wrong distance.\n");
    size_t k(0), l(0);
    for(size_t i(0); i < 3; ++i) 
      if(eq(std::abs(i_first->pos(i)), 0.5)) ++k;
      else if(is_null(i_first->pos(i))) ++l;
    LADA_DOASSERT(k == 2 and l == 1, "Wrong position.\n");
  }
  for(size_t i(0); i < 12; ++i, ++i_first)
  {
    LADA_DOASSERT(i_first != i_end, "Premature death.\n");
    LADA_DOASSERT(eq(i_first->distance, sqrt(0.75*0.75+2*0.25*0.25)), "Wrong distance.\n");
    size_t k(0), l(0);
    for(size_t i(0); i < 3; ++i) 
      if(eq(std::abs(i_first->pos(i)), 0.25)) ++k;
      else if(eq(std::abs(i_first->pos(i)), 0.75)) ++l;
    LADA_DOASSERT(k == 2 and l == 1, "Wrong position.\n");
  }
  for(size_t i(0); i < 6; ++i, ++i_first)
  {
    LADA_DOASSERT(i_first != i_end, "Premature death.\n");
    LADA_DOASSERT(eq(i_first->distance, 1e0), "Wrong distance.\n");
    size_t k(0), l(0);
    for(size_t i(0); i < 3; ++i) 
      if(eq(std::abs(i_first->pos(i)), 0e0)) ++k;
      else if(eq(std::abs(i_first->pos(i)), 1e0)) ++l;
    LADA_DOASSERT(k == 2 and l == 1, "Wrong position.\n");
  }
  for(size_t i(0); i < 12; ++i, ++i_first)
  {
    LADA_DOASSERT(i_first != i_end, "Premature death.\n");
    LADA_DOASSERT(eq(i_first->distance, sqrt(0.25*0.25+2*0.75*0.75)), "Wrong distance.\n");
    size_t k(0), l(0);
    for(size_t i(0); i < 3; ++i) 
      if(eq(std::abs(i_first->pos(i)), 0.75)) ++k;
      else if(eq(std::abs(i_first->pos(i)), 0.25)) ++l;
    LADA_DOASSERT(k == 2 and l == 1, "Wrong position.\n");
  }
  LADA_DOASSERT(i_first == i_end, "Unexpected life.\n");
}
int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  using namespace LaDa::math;
  typedef Structure< LADA_TYPE > t_Str; 
  Structure< LADA_TYPE > lattice;
  Neighbors neighbors(36);
  lattice.set_cell(0,0.5,0.5)
                  (0.5,0,0.5)
                  (0.5,0.5,0);
  lattice.add_atom(0,0,0, "Si")
                  (0.25,0.25,0.25, "Ge");
  typename t_Str::const_iterator i_atom = lattice.begin();
  typename t_Str::const_iterator i_atom_end = lattice.end();
  for(; i_atom != i_atom_end; ++i_atom)
  { 
    neighbors.origin = i_atom->pos();
    Neighbors::const_iterator i_first = neighbors.begin(lattice);
    Neighbors::const_iterator i_end = neighbors.end();
    check(i_first, i_end);
  }

  math::rMatrix3d cell;
  cell << 1, 1, 0, -5, 2, 0, 0, 0, 1;
  lattice = supercell(lattice, lattice.cell() * cell);
  i_atom = lattice.begin();
  i_atom_end = lattice.end();
  for(; i_atom != i_atom_end; ++i_atom)
  { 
    neighbors.origin = i_atom->pos();
    Neighbors::const_iterator i_first = neighbors.begin(lattice);
    Neighbors::const_iterator i_end = neighbors.end();
    check(i_first, i_end);
  }

  return 0;
}
