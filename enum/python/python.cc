//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/module.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/borrowed.hpp>

#include "find_all_cells.hpp"
#include "label_exchange.hpp"
#include "translation.hpp"
#include "transform.hpp"
#include "bitset.hpp"

BOOST_PYTHON_MODULE(enumeration)
{
  // loads lada.math first
  namespace bp = boost::python;
  bp::handle<> math( bp::borrowed(PyImport_ImportModule("lada.math")) );

  LaDa::Python::expose_find_all_cells();
  LaDa::Python::expose_label_exchange();
  LaDa::Python::expose_translation();
  LaDa::Python::expose_transform();
  LaDa::Python::expose_bitset();
}
