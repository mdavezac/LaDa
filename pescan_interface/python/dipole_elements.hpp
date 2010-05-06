#ifndef LADA_PYTHON_DIPOLE_ELEMENTS_HPP
#define LADA_PYTHON_DIPOLE_ELEMENTS_HPP

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LaDa
{
  namespace python
  {
    void expose_dipoles();
    void expose_oscillator_strength();
  }
} // namespace LaDa
#endif
