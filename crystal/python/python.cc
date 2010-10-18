#include "LaDaConfig.h"

#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/borrowed.hpp>
#include <boost/python/def.hpp>

#include <math/eigen.h>
#include <opt/types.h>

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

LaDa::types::t_real third_order(LaDa::math::rMatrix3d const & _matrix, LaDa::types::t_int _n)
{
  typedef LaDa::math::rVector3d rVector3d;
  typedef LaDa::types::t_real t_real;
  t_real result = 0e0;
  t_real ninv = 1e0 / t_real(_n);
  rVector3d d(-0.5,0,0);
  for(size_t i(0); i < _n; ++i, d(0) += ninv)
  {
    d(1) = -0.5;
    for(size_t j(0); j < _n; ++j, d(1) += ninv)
    {
      d(2) = -0.5;
      for(size_t k(0); k < _n; ++k, d(2) += ninv)
      {
        t_real min_dist = d.squaredNorm();
        for(int l(-1); l < 2; ++l)
          for(int m(-1); m < 2; ++m)
            for(int n(-1); n < 2; ++n)
            {
              t_real const m = (_matrix * (rVector3d(l, m, n) + d)).squaredNorm();
              if( m < min_dist ) min_dist = m;
            }
        result += min_dist;
      }
    }
  }
  return result  / (_matrix.determinant() * t_real(_n*_n*_n));
}

BOOST_PYTHON_MODULE(_crystal)
{
  namespace bp = boost::python;
  bp::scope scope;
  scope.attr("__doc__") = "This namespace is imported into lada.crystal.\n";
  bp::docstring_options doc_options(true, false);

  // loads lada.math first
  namespace bp = boost::python;
  bp::handle<> math( bp::borrowed(PyImport_ImportModule("lada.math")) );

  bp::def("third_order", &third_order);

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
