#ifndef LADA_MATH_GRUBER_H
#define LADA_MATH_GRUBER_H

#include "LaDaConfig.h"

#include "eigen.h"

namespace LaDa
{
  namespace math
  {
    //! Computes the Niggli cell of a lattice.
    rMatrix3d gruber(rMatrix3d const &_in, size_t itermax = 0);
    namespace details
    {
      //! function which needs to be compiled without optimization.
      bool no_opt_change_test(types::t_real _new, types::t_real _last);
    }
  }
}

#endif

