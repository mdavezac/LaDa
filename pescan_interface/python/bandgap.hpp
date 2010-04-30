//
//  Version: $Id$
//

#ifndef _LADA_PYTHON_BANDGAP_HPP_
#define _LADA_PYTHON_BANDGAP_HPP_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LaDa
{
  namespace python
  {
    void expose_bands();
    void expose_bandgap();
    void expose_oscillator_strength();
  }
} // namespace LaDa
#endif
