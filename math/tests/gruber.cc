#include "PyladaConfig.h"

#include "../math.h"

PyObject* testme(PyObject *_module, PyObject *)
{
  using namespace std;
  using namespace Pylada;
  using namespace Pylada::math;

  rMatrix3d cell, mult;
  cell << 0, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0;
  types::t_int const lim(PYLADA_LIM);
  Eigen::Matrix<types::t_real, 6, 1> params;
  params << 0.5, 0.5, 0.5, 0.5, 0.5, 0.5;

  if(not eq(cell, gruber(cell)))
  {
    PYLADA_PYERROR(internal, "Did not leave cell as is.\n");
    return NULL;
  }
  if(not eq(gruber_parameters(cell), params) )
  {
    PYLADA_PYERROR(internal, "Did not leave cell as is.\n");
    return NULL;
  }

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
    if(not is_integer(gruberred * cell.inverse()))
    {
      PYLADA_PYERROR(internal, "not a sublattice.\n");
      return NULL;
    }
    if(not is_integer(cell * gruberred.inverse()))
    {
      PYLADA_PYERROR(internal, "not a sublattice.\n");
      return NULL;
    }
    if(not eq(gruber_parameters(cell), params))
    {
      PYLADA_PYERROR(internal, "not a sublattice.\n");
      return NULL;
    }
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

PyMODINIT_FUNC init_gruber(void) 
{
  PyObject* module = Py_InitModule("_gruber", methods);
  if(not module) return;
  if(not Pylada::math::import()) return;
}
