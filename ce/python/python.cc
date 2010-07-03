#include "LaDaConfig.h"

#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/borrowed.hpp>

#include "ce.hpp"
#include "clusters.hpp"
#include "create_pairs.hpp"
#include "create_clusters.hpp"
#include "find_pis.hpp"
#include "mlcluster.hpp"
#include "mlclusters.hpp"

BOOST_PYTHON_MODULE(_ce)
{
  // loads lada.math first
  namespace bp = boost::python;
  bp::scope scope;
  scope.attr("__doc__") = "This namespace is imported into lada.crystal.\n";
  bp::docstring_options doc_options(true, false);
  bp::handle<> math( bp::borrowed(PyImport_ImportModule("lada.math")) );

  LaDa::Python::expose_ce();
  LaDa::Python::expose_clusters();
  LaDa::Python::expose_create_pairs();
  LaDa::Python::expose_create_clusters();
  LaDa::Python::expose_find_pis();
  LaDa::Python::expose_mlcluster();
  LaDa::Python::expose_mlclusters();
}
