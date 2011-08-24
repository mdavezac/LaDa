#include "LaDaConfig.h"
#include<iostream>
#include<string>

#include <opt/debug.h>
#include "../gruber.h"
#include "../misc.h"

using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::math;

  rMatrix3d cell, mult;
  cell << 0, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0;
  types::t_int const lim(LADA_LIM);
  Eigen::Matrix<types::t_real, 6, 1> params;
  params << 0.5, 0.5, 0.5, 0.5, 0.5, 0.5;

  LADA_DOASSERT(eq(cell, gruber(cell)), "Did not leave cell as is.\n");
  LADA_DOASSERT(eq(gruber_parameters(cell), params), "Did not leave cell as is.\n");

  bool cont = true;
  for(types::t_int a00(-lim); a00 <= lim and cont; ++a00)
  for(types::t_int a01(-lim); a01 <= lim and cont; ++a01)
  for(types::t_int a02(-lim); a02 <= lim and cont; ++a02)
  for(types::t_int a10(-lim); a10 <= lim and cont; ++a10)
  for(types::t_int a11(-lim); a11 <= lim and cont; ++a11)
  for(types::t_int a12(-lim); a12 <= lim and cont; ++a12)
  for(types::t_int a20(-lim); a20 <= lim and cont; ++a20)
  for(types::t_int a21(-lim); a21 <= lim and cont; ++a21)
  for(types::t_int a22(-lim); a22 <= lim and cont; ++a22)
  {
    mult << a00, a01, a02, a10, a11, a12, a20, a21, a22;
    if( neq(mult.determinant(), 1e0) ) continue;
    rMatrix3d const gruberred = gruber(cell*mult);
    LADA_DOASSERT(is_integer(gruberred * cell.inverse()), "not a sublattice.\n")
    LADA_DOASSERT(is_integer(cell * gruberred.inverse()), "not a sublattice.\n")
    LADA_DOASSERT(eq(gruber_parameters(cell), params), "not a sublattice.\n")
  }
  
  return 0;
}
