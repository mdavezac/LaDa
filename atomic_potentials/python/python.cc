//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/module.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/borrowed.hpp>

#include "bases.hpp"
#include "representation.hpp"
#include "sum_of_separables.hpp"
#include "separable.hpp"
#include "functions.hpp"
#include "collapse.hpp"
#include "fitting_set.hpp"
#include "values.hpp"

BOOST_PYTHON_MODULE(potentials)
{
  // loads lada.math first
  namespace bp = boost::python;
  bp::handle<> math( bp::borrowed(PyImport_ImportModule("lada.math")) );

  LaDa::Python::expose_representation();
  LaDa::Python::expose_bases();
  LaDa::Python::expose_sumofseps();
  LaDa::Python::expose_separable();
  LaDa::Python::expose_functions();
  LaDa::Python::expose_fittingset();
  LaDa::Python::expose_values();
  LaDa::Python::expose_collapse();
}
