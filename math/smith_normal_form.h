#ifndef _LADA_OPT_SMITH_NORMAL_FORM_H_
#define _LADA_OPT_SMITH_NORMAL_FORM_H_

#include "LaDaConfig.h"

#include "eigen.h"

namespace LaDa 
{
  namespace math
  {
    //! Computes smith normal form of a matrix \a  _S = \a _L \a _M \a _R.
    void smith_normal_form( math::iMatrix3d& _S, math::iMatrix3d & _L,
                            const math::iMatrix3d& _M, math::iMatrix3d &_R );
  } // namespace opt

} // namespace LaDa

#endif
