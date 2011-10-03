#include "LaDaConfig.h"

#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/borrowed.hpp>
#include <boost/python/def.hpp>

#include <numpy/arrayobject.h>

// #include <math/eigen.h>
// #include <opt/types.h>
// #include <physics/physics.h>
#include <python/std_vector.hpp>

#include "atom.hpp"
#include "set.hpp"
#include "structure.hpp"
// #include "read_structure.hpp"
// #include "enumerate.hpp"
// #include "smith.hpp"
// #include "layerdepth.hpp"
// #include "neighbors.hpp"
// #include "symmetry_operator.hpp"
// #include "which_site.hpp"
// #include "misc.hpp"
// #include "periodic_dnc.hpp"

// LaDa::types::t_real third_order(LaDa::math::rMatrix3d const & _matrix, LaDa::types::t_int _n)
// {
//   typedef LaDa::math::rVector3d rVector3d;
//   typedef LaDa::types::t_real t_real;
//   t_real result = 0e0;
//   t_real const ninv = 1e0 / t_real(_n);
//   t_real const maxdist = (_matrix * 10e0 * rVector3d(1,1,1)).squaredNorm();
//   for(size_t i(0); i < _n; ++i)
//   {
//     for(size_t j(0); j < _n; ++j)
//     {
//       for(size_t k(0); k < _n; ++k)
//       {
//         t_real min_dist = maxdist;
//         for(int l(-1); l < 2; ++l)
//           for(int m(-1); m < 2; ++m)
//             for(int n(-1); n < 2; ++n)
//             {
//               t_real const q = (_matrix * rVector3d( i*ninv+l-0.5, 
//                                                      j*ninv+m-0.5, 
//                                                      k*ninv+n-0.5  ) ).squaredNorm();
//               if( q < min_dist ) min_dist = q;
//             }
//         result += min_dist;
//       }
//     }
//   }
//   return result / (_matrix.determinant() * t_real(_n*_n*_n));
// }
//
// LaDa::types::t_unsigned nb_valence_states( LaDa::Crystal::TStructure<std::string> const &_str ) 
// {
//   LaDa::Crystal::TStructure<std::string>::t_Atoms::const_iterator i_atom = _str.atoms.begin();
//   LaDa::Crystal::TStructure<std::string>::t_Atoms::const_iterator i_atom_end = _str.atoms.end();
//   LaDa::types::t_unsigned bgstates = 0;
//   for(; i_atom != i_atom_end; ++i_atom)
//     bgstates += LaDa::Physics::Atomic::Charge( i_atom->type );
//   return bgstates;
// }

BOOST_PYTHON_MODULE(_crystal)
{
  namespace bp = boost::python;
  bp::scope scope;
  scope.attr("__docformat__") = "restructuredtext en";
  scope.attr("__doc__") = "This namespace is imported into lada.crystal.\n";
  bp::docstring_options doc_options(true, false);

  // imports numpy and loads lada.math first
  import_array();
  namespace bp = boost::python;
  bp::handle<> math( bp::borrowed(PyImport_ImportModule("lada.math")) );
  bp::handle<> error( bp::borrowed(PyImport_ImportModule("lada.error")) );

// bp::def("third_order", &third_order, (bp::arg("structure"), bp::arg("n")=200));
// bp::def( "nb_valence_states", &nb_valence_states, bp::arg("structure"), 
//          "Returns the number of `escan` valence states in a structure." );

  LaDa::python::expose_vector<std::string>("VectorStr", "Interface to Cpp vectors of string.");
  LaDa::python::expose_set();
  LaDa::python::expose_atom();
// LaDa::python::expose_structure();
// LaDa::Python::expose_structure();
// LaDa::Python::expose_lattice();
// LaDa::Python::expose_read_structure();
// LaDa::Python::expose_enumerate();
// LaDa::Python::expose_smith();
// LaDa::Python::expose_layerdepth();
// LaDa::Python::expose_neighbors();
// LaDa::Python::expose_symmetry_operator();
// LaDa::Python::expose_which_site();
// LaDa::python::expose_misc();
// LaDa::python::expose_periodic_dnc();
}
