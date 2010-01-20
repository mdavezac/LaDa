//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_value_policy.hpp>

#include "avogadro.hpp"
#include "matrix.hpp"
#include "../eigen.h"

BOOST_PYTHON_MODULE(math)
{
  namespace bp = boost::python;
  bp::scope scope;
  scope.attr("__doc__") = "This namespace defines mathematical types\n\nIt "
   "imports automatic conversions between numpy arrays and 3d-vectors and "
   "matrices.";
  bp::docstring_options doc_options(true, false);

  LaDa::Python::expose_eigen_vectors();
  LaDa::Python::expose_eigen_matrices();
}
