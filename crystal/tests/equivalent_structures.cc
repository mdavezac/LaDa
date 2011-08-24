#include "LaDaConfig.h"
#include<iostream>
#include<string>
#include<vector>

#include "../supercell.h"
#include "../equivalent_structures.h"

using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  using namespace LaDa::math;
  TemplateStructure< std::vector<std::string> > A; 
  A.set_cell(0,0.5,0.5)
            (0.5,0,0.5)
            (0.5,0.5,0);
  A.add_atom(0,0,0, "Si")
   .add_atom(0.25,0.25,0.25, "Si", "Ge");
  TemplateStructure< std::vector<std::string> > B = A.copy();

  LADA_DOASSERT(equivalent(A, A, true, 1e-5), "A is not equivalent to itself.\n");
  LADA_DOASSERT(equivalent(A, B, true, 1e-5), "A is not equivalent to itself.\n");
  B.scale() = 3.0;
  LADA_DOASSERT(not equivalent(A, B, true, 1e-5), "A is equivalent to B.\n");
  LADA_DOASSERT(equivalent(A, B, false, 1e-5), "A is not equivalent to B.\n");
  B.cell() *= 0.5;
  foreach(TemplateStructure< std::vector<string> >::reference b, B) b.pos *= 0.5;
  LADA_DOASSERT(not equivalent(A, B, true, 1e-5), "A is equivalent to B.\n");
  LADA_DOASSERT(equivalent(A, B, false, 1e-5), "A is not equivalent to B.\n");

  B = A.copy();
  B[1].pos *= -1;
  LADA_DOASSERT(equivalent(A, B, true, 1e-5), "A is not equivalent to B.\n");

  math::Affine3d affine(math::AngleAxis(0.5, math::rVector3d::UnitX()));
  B = A.transform(affine);
  LADA_DOASSERT(equivalent(A, B, true, 1e-5), "A is not equivalent to B.\n");
  B[0].type.push_back("Ge");
  B[1].type.pop_back(); 
  LADA_DOASSERT(not equivalent(A, B, true, 1e-5), "A is equivalent to B.\n");
  affine = math::Translation(0.5, 0, 0);
  B = A.transform(affine);
  LADA_DOASSERT(equivalent(A, B, true, 1e-5), "A is not equivalent to B.\n");
  B.scale() = 3e0;
  B.cell() *= 0.5;
  foreach(TemplateStructure< std::vector<string> >::reference b, B) b.pos *= 0.5;
  LADA_DOASSERT(not equivalent(A, B, true, 1e-5), "A is equivalent to B.\n");
  LADA_DOASSERT(equivalent(A, B, false, 1e-5), "A is not equivalent to B.\n");


  return 0;
}
