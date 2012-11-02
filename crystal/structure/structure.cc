#include "LaDaConfig.h"

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL lada_math_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>

#include "structure.h"

namespace LaDa
{

  namespace crystal 
  {
    // Transforms a structure according to an affine transformation.
    void Structure::transform(Eigen::Matrix<types::t_real, 4, 3> const &_affine)
    {
      cell() = _affine.block<3,3>(0,0) * cell();
      iterator i_first = begin();
      iterator const i_end = end();
      for(; i_first != i_end;  ++i_first)
        i_first->pos() = _affine.block<3,3>(0, 0) * i_first->pos() + ~_affine.block<1, 3>(3, 0);
    }
  } // namespace Crystal

} // namespace LaDa
