#include "LaDaConfig.h"
#include<iostream>
#include<string>
#include<vector>

#include "../supercell.h"
#include "../equivalent_structures.h"

using namespace LaDa;
using namespace LaDa::crystal;
using namespace LaDa::math;

void test(TemplateStructure< std::vector<std::string> > const _A, 
          TemplateStructure< std::vector<std::string> > const _B)
{
  TemplateStructure< std::vector<std::string> > B = _B.copy();
  LADA_DOASSERT(equivalent(_A, B, true, 1e-5), "A is not equivalent to B.\n");
  B.scale() = 3.0;
  LADA_DOASSERT(not equivalent(_A, B, true, 1e-5), "A is equivalent to B.\n");
  LADA_DOASSERT(equivalent(_A, B, false, 1e-5), "A is not equivalent to B.\n");
  B.cell() *= 0.5;
  foreach(TemplateStructure< std::vector<std::string> >::reference b, B) b.pos *= 0.5;
  LADA_DOASSERT(not equivalent(_A, B, true, 1e-5), "A is equivalent to B.\n");
  LADA_DOASSERT(equivalent(_A, B, false, 1e-5), "A is not equivalent to B.\n");
}
int main()
{
  TemplateStructure< std::vector<std::string> > A, B; 
  A.set_cell(0,0.5,0.5)
            (0.5,0,0.5)
            (0.5,0.5,0);
  A.add_atom(0,0,0, "Si")
   .add_atom(0.25,0.25,0.25, "Si", "Ge");

  test(A, A);
 
  B = A.copy();
  B[1].pos *= -1;
  test(A, B);
 
  math::Affine3d affine(math::AngleAxis(0.5, math::rVector3d::UnitX()));
  B = A.transform(affine);
  test(A, B);
 
  B[0].type.push_back("Ge");
  B[1].type.pop_back(); 
  test(A, B);
 
  affine = math::Translation(0.5, 0, 0);
  B = A.transform(affine);
  test(A, B);
 
  affine = math::AngleAxis(0.5 * math::pi, math::rVector3d::UnitX());
  B = A.transform(affine);
  test(A, B);
 
  affine = math::AngleAxis(0.27 * math::pi, math::rVector3d::UnitX());
  B = A.transform(affine);
  test(A, B);
 
  affine = math::AngleAxis(0.27 * math::pi, math::rVector3d::UnitX()) * math::Translation(0.5,0,0);
  B = A.transform(affine);
  test(A, B);
 
  affine = math::AngleAxis(0.27 * math::pi, math::rVector3d::UnitX());
  affine = math::Translation(0.5,0,0) * affine;
  B = A.transform(affine);
  test(A, B);
 
  B = A.copy();
  B.set_cell(-0.5, -0.5, 0)
            (-0.5, 0, -0.5)
            (0, -0.5, -0.5);
  test(A, B);
  B[1].pos = rVector3d(-0.75, -0.75, -0.75);
  test(A, B);

  A[0].type.clear(); A[0].type.push_back("Si");
  A[1].type.clear(); A[1].type.push_back("Si");
  rMatrix3d cell;
  cell << 3, 0, 0, 0, 0.5, -0.5, 0, 0.5, 0.5;
  A = supercell(A, cell);
  A[0].type[0] = "Ge";
  B = A.copy();
  B[0].type[0] = "Si";
  B[5].type[0] = "Ge";
  test(A, B);

  A[0].type[0] = "Ge";
  A[2].type[0] = "Ge";
  A[3].type[0] = "Ge";

  B[5].type[0] = "Ge";
  B[6].type[0] = "Si";
  B[7].type[0] = "Ge";
  B[8].type[0] = "Ge";
  test(A, B);



  return 0;
}
