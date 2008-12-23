//
//  Version: $Id$
//

#ifndef _LADA_PYTHON_VFF_HPP_
#define _LADA_PYTHON_VFF_HPP_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <crystal/structure.h>

namespace LaDa
{
  namespace Python
  {
    void expose_vff();
    void expose_layered_vff();
  }
} // namespace LaDa
#endif
