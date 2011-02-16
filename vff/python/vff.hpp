#ifndef _LADA_PYTHON_VFF_HPP_
#define _LADA_PYTHON_VFF_HPP_

#include "LaDaConfig.h"

#include <crystal/structure.h>

namespace LaDa
{
  namespace python
  {
    void expose_vff();
    void expose_layeredvff();
  }
} // namespace LaDa
#endif
