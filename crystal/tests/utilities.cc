#include "LaDaConfig.h"
#include<iostream>
#include<string>
#include<vector>

#include "../utilities.h"

using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  using namespace LaDa::math;
  rMatrix3d matrix;
  matrix << 0, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0;
  rVector3d vec(0.75, -0.25, 0);
  LADA_ASSERT(eq(into_cell(vec, matrix), rVector3d(0.25, 0.25, 0)), "Not in cell.\n")
  LADA_ASSERT(eq(into_voronoi(vec, matrix), rVector3d(0.25, 0.25, 0)), "Not in cell.\n")

  matrix << 0, 0.5, 1.5,
            0.5, 0, 1.5, 
            0.5, 0.5, 0;
  vec = matrix*vec;
  LADA_ASSERT(eq(into_cell(vec, matrix), rVector3d(0.375, 0.375, 0.75)), "Not in cell.\n")
  LADA_ASSERT(eq(into_voronoi(vec, matrix), rVector3d(-0.125, -0.125, -0.25)), "Not in cell.\n")
  LADA_ASSERT(is_integer( matrix.inverse() * 
                            (into_voronoi(vec, matrix) - into_cell(vec, matrix)) ), "Not in cell.\n")
  return 0;
}
