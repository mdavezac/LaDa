//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>
#include <boost/python/module.hpp>

#include "clj.hpp"

BOOST_PYTHON_MODULE(models)
{
  LaDa::Python::expose_clj();
}
