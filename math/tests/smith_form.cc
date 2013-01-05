#include "PyladaConfig.h"

#include<iostream>

#include <cstdlib>
#include <time.h>

#include "../math.h"

#define PYLADA_DOASSERT(a,b)                \
        {                                 \
          if((not (a)))                   \
          {                               \
            PYLADA_PYERROR(internal, b);    \
            return NULL;                  \
          }                               \
        }

PyObject* testme(PyObject* _module, PyObject *)
{
  using namespace std;
  using namespace Pylada;
  using namespace Pylada::math;

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
    for(int j(0); j < cell.rows(); ++j)
      for(int k(0); k < cell.cols(); ++k)
        if(k != j) { PYLADA_DOASSERT(smith(j,k) == 0, "Non-zero off diagonal.\n") }
        else { PYLADA_DOASSERT(smith(j,k) != 0, "Zero on diagonal.\n") }
    for(int j(0); j < cell.rows() - 1; ++j)
      PYLADA_DOASSERT(smith(j+1,j+1) % smith(j,j) == 0, "Not a factor.\n")
    PYLADA_DOASSERT( smith == left * cell * right, "Not a transform.\n");
    PYLADA_DOASSERT( std::abs(left.determinant()) > 1e-12, "Left matrix not invertible.\n")
    PYLADA_DOASSERT( std::abs(right.determinant()) > 1e-12, "Right matrix not invertible.\n")
  }
  Py_RETURN_TRUE;
}

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

#ifdef PYLADA_DECLARE
#  error PYLADA_DECLARE already defined.
#endif
#define PYLADA_DECLARE(name, args) {#name, (PyCFunction)name, METH_ ## args, ""} 

static PyMethodDef methods[] = { 
  PYLADA_DECLARE(testme, NOARGS),
  {NULL},
};

#undef PYLADA_DECLARE

PyMODINIT_FUNC init_smith(void) 
{
  PyObject* module = Py_InitModule("_smith", methods);
  if(not module) return;
  if(not Pylada::math::import()) return;
}
