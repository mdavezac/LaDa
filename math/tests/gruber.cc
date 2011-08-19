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
  LADA_ASSERT(eq(cell, gruber(cell)), "Did not leave cell as is.\n");
  types::t_int const lim(2);
  for(types::t_int a00(-lim); a00 <= lim; ++a00)
  for(types::t_int a01(-lim); a01 <= lim; ++a01)
  for(types::t_int a02(-lim); a02 <= lim; ++a02)
  for(types::t_int a10(-lim); a10 <= lim; ++a10)
  for(types::t_int a11(-lim); a11 <= lim; ++a11)
  for(types::t_int a12(-lim); a12 <= lim; ++a12)
  for(types::t_int a20(-lim); a20 <= lim; ++a20)
  for(types::t_int a21(-lim); a21 <= lim; ++a21)
  for(types::t_int a22(-lim); a22 <= lim; ++a22)
  {
    mult << a00, a01, a02, a10, a11, a12, a20, a21, a22;
//   if( not is_identity(mult * (~mult)) ) continue;
    if( neq(mult.determinant(), 1e0) ) continue;
    std::cout << gruber(mult) << "\n\n";
  }
  
  std::cout << gruber(cell) << "\n";
  
  return 0;
}
