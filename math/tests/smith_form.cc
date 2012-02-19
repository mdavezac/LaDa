#include "LaDaConfig.h"

#include<iostream>

#include <opt/debug.h>
#include <cstdlib>
#include <time.h>

#include "../smith_normal_form.h"
#include "../fuzzy.h"


using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::math;

  types::t_int const seed = time(NULL); // 1318126056; //
  std::cout << seed << "\n";
  srand(seed);

  size_t const n = 5;
  for(size_t i(0); i < 10000; ++i)
  {
    iMatrix3d cell, left, right, smith;
    do
    {
      cell << rand()%(n<<1)-n, rand()%(n<<1)-n, rand()%(n<<1)-n,
              rand()%(n<<1)-n, rand()%(n<<1)-n, rand()%(n<<1)-n,
              rand()%(n<<1)-n, rand()%(n<<1)-n, rand()%(n<<1)-n;
    } while( is_null(cell.determinant()) );
    smith_normal_form(smith, left, cell, right);
    for(size_t j(0); j < cell.rows(); ++j)
      for(size_t k(0); k < cell.cols(); ++k)
        if(k != j) { LADA_DOASSERT(smith(j,k) == 0, "Non-zero off diagonal.\n") }
        else { LADA_DOASSERT(smith(j,k) != 0, "Zero on diagonal.\n") }
    for(size_t j(0); j < cell.rows() - 1; ++j)
      LADA_DOASSERT(smith(j+1,j+1) % smith(j,j) == 0, "Not a factor.\n")
    LADA_DOASSERT( smith == left * cell * right, "Not a transform.\n");
    LADA_DOASSERT( std::abs(left.determinant()) > 1e-12, "Left matrix not invertible.\n")
    LADA_DOASSERT( std::abs(right.determinant()) > 1e-12, "Right matrix not invertible.\n")
  }
  return 0;
}
