#include "LaDaConfig.h"


#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <time.h>

#include <math/smith_normal_form.h>

#include "../supercell.h"
#include "../space_group.h"
#include "../equivalent_structures.h"

using namespace LaDa;
using namespace LaDa::crystal;
using namespace LaDa::math;
typedef TemplateStructure< std::vector<std::string> > t_Str;
#define LADA_RAND(s) types::t_real(rand()%s)/types::t_real(s)

void scale(t_Str const &_A, t_Str const &_B)
{
  t_Str B = _B.copy();
  LADA_DOASSERT(equivalent(_A, B, true, 1e-5), "A is not equivalent to B.\n");
  B.scale() = 3.0;
  LADA_DOASSERT(not equivalent(_A, B, true, 1e-5), "A is equivalent to B.\n");
  LADA_DOASSERT(equivalent(_A, B, false, 1e-5), "A is not equivalent to B.\n");
  B.cell() *= 0.5;
  foreach(t_Str::reference b, B) b.pos *= 0.5;
  LADA_DOASSERT(not equivalent(_A, B, true, 1e-5), "A is equivalent to B.\n");
  LADA_DOASSERT(equivalent(_A, B, false, 1e-5), "A is not equivalent to B.\n");
}
void motif(t_Str const &_A, t_Str const &_B)
{
  boost::shared_ptr<t_SpaceGroup> sg = cell_invariants(_A.cell());
  t_SpaceGroup::const_iterator i_first = sg->begin();
  t_SpaceGroup::const_iterator const i_end = sg->begin();
  for(; i_first != i_end; ++i_first)
  {
    t_Str B = _B.copy();
    foreach(t_Str::reference b, B)
    {
      b.pos = i_first->linear() * b.pos;
      scale(_A, B);
    }
  }
}
void basis(t_Str const &_A, t_Str const &_B)
{
  math::Affine3d affine(math::AngleAxis(0.5*math::pi, math::rVector3d::UnitX()));
  t_Str B = _B.transform(affine);
  motif(_A, B);
  affine = math::AngleAxis(-math::pi, math::rVector3d::UnitX());
  B = _B.transform(affine);
  motif(_A, B);
  affine = math::AngleAxis(-0.13*math::pi, math::rVector3d::UnitX());
  B = _B.transform(affine);
  motif(_A, B);
  affine = math::Translation(0.25, 0.25, 0.25);
  B = _B.transform(affine);
  motif(_A, B);

  affine =   math::AngleAxis(LADA_RAND(100)*2.0*math::pi, math::rVector3d::UnitX())
           * math::Translation(LADA_RAND(100)-0.5, LADA_RAND(100)-0.5, LADA_RAND(100)-0.5);
  B = _B.transform(affine);
  motif(_A, B);
}
void decoration(t_Str const &_A, t_Str const &_B, t_Str const _latt)
{
  // map A's atoms with linear smith index.
  t_SmithTransform st = smith_transform(_latt.cell(), _A.cell());
  std::vector<size_t> indices(_A.size());
  t_Str::const_iterator i_atom = _A.begin();
  t_Str::const_iterator const i_end = _A.end();
  for(size_t i(0); i_atom != i_end; ++i_atom, ++i)
    indices[linear_smith_index(st, i_atom->site, i_atom->pos - _latt[i_atom->site].pos)] = i;

  // loop over decoration translations
  for(i_atom = _A.begin(); i_atom != i_end; ++i_atom)
  {
    if(i_atom->site != _A[0].site) continue; // only primitive lattice translations.
    math::rVector3d const trans = i_atom->pos - _A[0].pos;
    t_Str B = _A.copy();
    std::vector<size_t>::const_iterator i_ind = indices.begin();
    std::vector<size_t>::const_iterator const i_ind_end = indices.end();
    for(; i_ind != i_ind_end; ++i_ind)
    {
      math::rVector3d const vec = _A[*i_ind].pos + trans - _latt[_A[*i_ind].site].pos;
      if(not is_integer(_latt.cell().inverse() * vec))
      {
        std::cout << ~(_latt.cell().inverse() * vec) << "\n";
        LADA_DOASSERT(false, "Not on lattice.\n" )
      }
      B[ indices[linear_smith_index(st, _A[*i_ind].site, vec)] ] = _A[*i_ind];
    }
    basis(_A, B);
  }
}

int main()
{
  types::t_int const seed = time(NULL);
  std::cout << seed << "\n";
  srand(seed);

  t_Str A;
  A.set_cell(0,0.5,0.5)
            (0.5,0,0.5)
            (0.5,0.5,0);
  A.add_atom(0,0,0, "Si")
            (0.25,0.25,0.25, "Si", "Ge");
  basis(A, A);

  A[0].type.clear(); A[0].type.push_back("Si");
  A[1].type.clear(); A[1].type.push_back("Si");
  t_Str latt = A.copy();
  rMatrix3d cell;
  cell << 3, 0, 0, 0, 0.5, -0.5, 0, 0.5, 0.5;
  A = supercell(latt, cell);
  A[0].type[0] = "Ge";
  A[1].type[0] = "Ge";
  A[3].type[0] = "Ge";

  decoration(A, A, latt);

  cell << 2, 0, 0, 0, 2, 0, 0, 0, 2;
  A = supercell(latt, cell);
  for(size_t i(0); i < 5; ++i)
  {
    t_Str::iterator i_atom = A.begin();
    t_Str::iterator const i_end = A.end();
    types::t_real const x = LADA_RAND(100) * 0.5;
    for(; i_atom != i_end; ++i_atom) i_atom->type[0] = LADA_RAND(100) > x ? "Si": "Ge";
    decoration(A, A, latt);
  }
 
  for(size_t i(0); i < 10; ++i)
  {
    do
    {
      cell << rand()%10-5, rand()%10-5, rand()%10-5,
              rand()%10-5, rand()%10-5, rand()%10-5,
              rand()%10-5, rand()%10-5, rand()%10-5;
    } while( is_null(cell.determinant()) );
    A = supercell(latt, latt.cell() * cell);
    t_Str::iterator i_atom = A.begin();
    t_Str::iterator const i_end = A.end();
    types::t_real const x = LADA_RAND(100) * 0.5;
    for(; i_atom != i_end; ++i_atom) i_atom->type[0] = LADA_RAND(100) > x ? "Si": "Ge";
    decoration(A, A, latt);
  }

  return 0;
}
