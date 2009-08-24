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

BOOST_PYTHON_MODULE(potentials)
{
  LaDa::Python::expose_representation();
  LaDa::Python::expose_bases();
  LaDa::Python::expose_collapse();
}
