//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>

#include "cgs.hpp"

BOOST_PYTHON_MODULE(Minimizer)
{
  LaDa::Python::expose_cgs();
  LaDa::Python::expose_llsq();
  LaDa::Python::expose_mul_mat_vec();
}
