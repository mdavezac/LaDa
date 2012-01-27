#ifndef LADA_OPT_SMITH_NORMAL_FORM_H
#define LADA_OPT_SMITH_NORMAL_FORM_H

#include "LaDaConfig.h"

#include "eigen.h"

namespace LaDa 
{
  namespace math
  {
    //! Computes smith normal form of a matrix \a  _S = \a _L \a _M \a _R.
    void smith_normal_form( iMatrix3d& _S, iMatrix3d & _L, 
                            const iMatrix3d& _M, iMatrix3d &_R );
  } // namespace opt
} // namespace LaDa

#endif
