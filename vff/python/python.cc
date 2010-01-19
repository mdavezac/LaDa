//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/module.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/borrowed.hpp>

#include "vff.hpp"

BOOST_PYTHON_MODULE(vff)
{
  // loads lada.math first
  namespace bp = boost::python;
  bp::handle<> math( bp::borrowed(PyImport_ImportModule("lada.math")) );

  LaDa::Python::expose_vff();
  LaDa::Python::expose_layeredvff();
}
