//
//  Version: $Id$
//

#ifndef _LADA_PYTHON_CGS_HPP_
#define _LADA_PYTHON_CGS_HPP_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LaDa
{
  namespace Python
  {
    void expose_cgs();
    void expose_llsq();
    void expose_mul_mat_vec();
  }
} // namespace LaDa
#endif
