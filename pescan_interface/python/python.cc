//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>

#include "escan.hpp"
#include "bandgap.hpp"
#include "emass.hpp"
#include "dipole_elements.hpp"

BOOST_PYTHON_MODULE(_escan)
{
  namespace bp = boost::python;
  bp::docstring_options doc_options(true, false);
  bp::scope scope;
  scope.attr("__doc__") = "This namespace is imported into lada.escan.\n";

  LaDa::Python::expose_escan_parameters();
  LaDa::Python::expose_escan();
  LaDa::Python::expose_bands();
  LaDa::Python::expose_bandgap();
  LaDa::Python::expose_oscillator_strength();
  LaDa::Python::expose_genpot();
  LaDa::Python::expose_emass();
  LaDa::Python::expose_dipoles();
}
