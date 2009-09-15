//
//  Version: $Id$
//

#ifndef _LADA_ENUM_PYTHON_OPERATIONS_HPP_
#define _LADA_ENUM_PYTHON_OPERATIONS_HPP_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LaDa
{
  namespace Python
  {
    void expose_label_exchange();
    void expose_translation();
    void expose_transform();
  }
} // namespace LaDa
#endif
