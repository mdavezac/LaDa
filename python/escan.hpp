//
//  Version: $Id$
//

#ifndef _LADA_PYTHON_ESCAN_HPP_
#define _LADA_PYTHON_ESCAN_HPP_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LaDa
{
  namespace Python
  {
    void expose_escan_parameters();
    void expose_escan();
  }
} // namespace LaDa
#endif
