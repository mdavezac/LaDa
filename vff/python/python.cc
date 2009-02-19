//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>

#include "vff.hpp"

BOOST_PYTHON_MODULE(Vff)
{
  LaDa::Python::expose_vff();
}
