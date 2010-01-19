//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/module.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/borrowed.hpp>

#include "lattice.hpp"
#include "structure.hpp"
#include "atom.hpp"
#include "read_structure.hpp"
#include "enumerate.hpp"
#include "smith.hpp"
#include "layerdepth.hpp"
#include "neighbors.hpp"
#include "symmetry_operator.hpp"
#include "which_site.hpp"

BOOST_PYTHON_MODULE(_crystal)
{
  // loads lada.math first
  namespace bp = boost::python;
  bp::handle<> math( bp::borrowed(PyImport_ImportModule("lada.math")) );

  LaDa::Python::expose_atom();
  LaDa::Python::expose_structure();
  LaDa::Python::expose_lattice();
  LaDa::Python::expose_read_structure();
  LaDa::Python::expose_enumerate();
  LaDa::Python::expose_smith();
  LaDa::Python::expose_layerdepth();
  LaDa::Python::expose_neighbors();
  LaDa::Python::expose_symmetry_operator();
  LaDa::Python::expose_which_site();
}
