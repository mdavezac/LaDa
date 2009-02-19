//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>

#include "escan.hpp"
#include "bandgap.hpp"

BOOST_PYTHON_MODULE(Pescan)
{
  LaDa::Python::expose_escan_parameters();
  LaDa::Python::expose_escan();
  LaDa::Python::expose_bands();
  LaDa::Python::expose_bandgap();
}
