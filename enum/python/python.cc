#include "LaDaConfig.h"

#include <boost/python/module.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/borrowed.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>
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
  bp::scope scope;
  scope.attr("__doc__") = "This namespace is imported into lada.crystal.\n";
  bp::docstring_options doc_options(true, false);
  bp::handle<> math( bp::borrowed(PyImport_ImportModule("lada.math")) );
  bp::handle<> crystal( bp::borrowed(PyImport_ImportModule("lada.crystal")) );

  LaDa::Python::expose_find_all_cells();
  LaDa::Python::expose_label_exchange();
  LaDa::Python::expose_translation();
  LaDa::Python::expose_transform();
  LaDa::Python::expose_bitset();
}
