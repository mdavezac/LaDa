//
//  Version: $Id$
//

#ifndef _LADA_OPT_SMITH_NORMAL_FORM_H_
#define _LADA_OPT_SMITH_NORMAL_FORM_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <atat/vectmac.h>
#include "types.h"

namespace LaDa 
{
  namespace opt
  {
    //! Computes smith normal form of a matrix \a  _S = \a _L \a _M \a _R.
    void smith_normal_form( atat::iMatrix3d& _S, atat::iMatrix3d & _L,
                            const atat::iMatrix3d& _M, atat::iMatrix3d &_R );
  } // namespace opt

} // namespace LaDa

#endif
