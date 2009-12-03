//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>

#include "escan.hpp"
#include "bandgap.hpp"
#include "emass.hpp"

BOOST_PYTHON_MODULE(_escan)
{
  LaDa::Python::expose_escan_parameters();
  LaDa::Python::expose_escan();
  LaDa::Python::expose_bands();
  LaDa::Python::expose_bandgap();
  LaDa::Python::expose_oscillator_strength();
  LaDa::Python::expose_genpot();
  LaDa::Python::expose_emass();
}
