//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>

#include "ce.hpp"

BOOST_PYTHON_MODULE(ce)
{
  LaDa::Python::expose_ce();
}