#include "LaDaConfig.h"
#include<iostream>
#include<string>
#include<vector>

#include "../supercell.h"
#include "../equivalent_structures.h"

using namespace LaDa;
using namespace LaDa::crystal;
using namespace LaDa::math;

void test(TemplateStructure< std::vector<std::string> > const _A)
{
  TemplateStructure< std::vector<std::string> > B = _A.copy();
  LADA_DOASSERT(equivalent(_A, B, true, true, 1e-5), "A is not equivalent to itself.\n");
  B.scale() = 3.0;
  LADA_DOASSERT(not equivalent(_A, B, true, true, 1e-5), "A is equivalent to B.\n");
  LADA_DOASSERT(equivalent(_A, B, false, true, 1e-5), "A is not equivalent to B.\n");
  B.cell() *= 0.5;
  foreach(TemplateStructure< std::vector<std::string> >::reference b, B) b.pos *= 0.5;
  LADA_DOASSERT(not equivalent(_A, B, true, true, 1e-5), "A is equivalent to B.\n");
  LADA_DOASSERT(equivalent(_A, B, false, true, 1e-5), "A is not equivalent to B.\n");
}
int main()
{
  TemplateStructure< std::vector<std::string> > A, B; 
  A.set_cell(0,0.5,0.5)
            (0.5,0,0.5)
            (0.5,0.5,0);
  A.add_atom(0,0,0, "Si")
   .add_atom(0.25,0.25,0.25, "Si", "Ge");

  test(A);

  B = A.copy();
  B[1].pos *= -1;
  test(B);

  math::Affine3d affine(math::AngleAxis(0.5, math::rVector3d::UnitX()));
  B = A.transform(affine);
  test(B);

  B[0].type.push_back("Ge");
  B[1].type.pop_back(); 
  LADA_DOASSERT(not equivalent(A, B, true, false, 1e-5), "A is equivalent to B.\n");
  test(B);

  affine = math::Translation(0.5, 0, 0);
  B = A.transform(affine);
  test(B);

  affine = math::AngleAxis(0.5 * math::pi, math::rVector3d::UnitX());
  B = A.transform(affine);
  test(B);

  affine = math::AngleAxis(0.27 * math::pi, math::rVector3d::UnitX());
  B = A.transform(affine);
  test(B);

  affine = math::AngleAxis(0.27 * math::pi, math::rVector3d::UnitX()) * math::Translation(0.5,0,0);
  B = A.transform(affine);
  test(B);

  affine = math::AngleAxis(0.27 * math::pi, math::rVector3d::UnitX());
  affine = math::Translation(0.5,0,0) * affine;
  B = A.transform(affine);
  test(B);

  return 0;
}
