//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>

#include "cgs.hpp"
#include "interpolated_gradient.hpp"
#include "minimizer.hpp"

BOOST_PYTHON_MODULE(minimizer)
{
  LaDa::Python::expose_cgs();
  LaDa::Python::expose_llsq();
  LaDa::Python::expose_interpolated_gradient();
  LaDa::Python::expose_minimizer();
}
