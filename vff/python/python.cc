//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>

#include "vff.hpp"

BOOST_PYTHON_MODULE(vff)
{
  LaDa::Python::expose_vff();
  LaDa::Python::expose_layeredvff();
}
