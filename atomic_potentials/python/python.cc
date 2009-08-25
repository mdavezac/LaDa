//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>

#include "bases.hpp"
#include "representation.hpp"
#include "sum_of_separables.hpp"
#include "collapse.hpp"
#include "fitting_set.hpp"
#include "values.hpp"

BOOST_PYTHON_MODULE(potentials)
{
  LaDa::Python::expose_representation();
  LaDa::Python::expose_bases();
// LaDa::Python::expose_sumofseps();
 // LaDa::Python::expose_fittingset();
 // LaDa::Python::expose_values();
 // LaDa::Python::expose_collapse();
}
