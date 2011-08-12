#include "LaDaConfig.h"

#include<iostream>

#include <opt/debug.h>

#include "../smith_normal_form.h"
#include "../fuzzy.h"
#include "../misc.h"


using namespace LaDa::math;
inline bool periodic_images( t_SmithTransform const &_transformation, 
                             rVector3d const &_a, rVector3d const &_b )
  { return eq(smith_index(_transformation, _a - _b), iVector3d::Zero()); }
using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::math;
  t_SmithTransform transform;
  rMatrix3d unitcell, supercell, testcell;
  unitcell << 0, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0;

  supercell << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  transform = smith_transform(unitcell, supercell);
  testcell << -1, 1, 1, 0, 2, 0, 0, 0, 2;
  LADA_ASSERT(eq(transform.get<0>(), testcell), "Incorrect transform.\n")
  LADA_ASSERT(eq(transform.get<1>(), iVector3d(1, 2, 2)), "Incorrect quotient.\n")
  LADA_ASSERT(eq(smith_index(transform, rVector3d(1, -0.5, 0.5)), iVector3d(0, 1, 1)),
               "Incorrect smith index.\n" )
  LADA_ASSERT(eq(smith_index(transform, rVector3d(2, 0.5, 0.5)), iVector3d(0, 1, 1)),
               "Incorrect smith index.\n" )
  LADA_ASSERT(periodic_images(transform, rVector3d(1, -0.5, 0.5), rVector3d(2, 0.5, 0.5)),
              "Are not perioric images.\n")
  LADA_ASSERT(periodic_images(transform, rVector3d(1, -0.425, 0.578), rVector3d(2, 0.575, 0.578)),
              "Are not perioric images.\n")

  bool result = false;
  try { periodic_images(transform, rVector3d(1, -0.575, 0.578), rVector3d(2, 0.575, 0.578)); }
  catch( error::off_lattice_position &_e ) { result = true; }
  LADA_ASSERT(result, "Did not catch exception.");
  result = true;
  try { periodic_images(transform, rVector3d(1, -0.575, 0.578), rVector3d(2, 0.575, 0.578)); }
  catch( error::unideal_lattice &_e ) { result = true; }
  LADA_ASSERT(result, "Did not catch exception.");

  supercell << 0, 0.4, 0.4, 0.5, 0, 0.5, 0.5, 0.5, 0;
  result = false;
  try { smith_transform(unitcell, supercell); }
  catch( error::not_a_supercell &_e ) { result = true; }
  catch( error::root &_e ) { LADA_ASSERT(true, "Should not be here.\n"); }
  LADA_ASSERT(result, "Did not catch exception.");
  result = false;
  try { smith_transform(unitcell, supercell); }
  catch( error::root &_e ) { result = true; }
  LADA_ASSERT(result, "Did not catch exception.");
  
  supercell << 3, 5, 0, 0, -1, 0, -2, -2, 1;
  supercell = unitcell * supercell;
  transform = smith_transform(unitcell, supercell);
  testcell << 0, 2, 0, 1, -1, 1, 5, 1, 3;
  LADA_ASSERT(eq(transform.get<0>(), testcell), "Incorrect transform.\n")
  LADA_ASSERT(eq(transform.get<1>(), iVector3d(1, 1, 3)), "Incorrect quotient.\n")
  LADA_ASSERT(eq(smith_index(transform, rVector3d(1, -0.5, 0.5)), iVector3d(0, 0, 0)),
               "Incorrect smith index.\n" )
  LADA_ASSERT(eq(smith_index(transform, rVector3d(-2, 0.5, 0.5)), iVector3d(0, 0, 1)),
               "Incorrect smith index.\n" )
  LADA_ASSERT(not periodic_images(transform, rVector3d(1, -0.5, 0.5), rVector3d(-2, 0.5, 0.5)),
              "Are not perioric images.\n")

  return 0;
}
