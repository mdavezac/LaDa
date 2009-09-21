//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <sstream>
#include <complex>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/errors.hpp>

#include <python/std_vector.hpp>
#include <python/misc.hpp>

#include "../structure.h"
#include "../lattice.h"
#include "../symmetry_operator.h"

#include "symmetry_operator.hpp"


namespace LaDa
{
  namespace Python
  {
    namespace bp = boost::python;

    boost::shared_ptr< std::vector<Crystal::SymmetryOperator> >
      get_point_group_symmetries( atat::rMatrix3d const &_cell, types::t_real _tolerance = -1e0 )
      {
        boost::shared_ptr< std::vector<Crystal::SymmetryOperator> > 
          result = Crystal::get_point_group_symmetries( _cell, _tolerance );
        if( not result )
        {
          PyErr_SetString(PyExc_RuntimeError, "Error while computing symmetry operations.\n");
          bp::throw_error_already_set();
        }
        return result;
      }
    boost::shared_ptr< std::vector<Crystal::SymmetryOperator> >
      get_space_group_symmetries( Crystal::Lattice const &_lattice, types::t_real _tolerance = -1e0 )
      {
        boost::shared_ptr< std::vector<Crystal::SymmetryOperator> > 
          result = Crystal::get_space_group_symmetries( _lattice, _tolerance );
        if( not result )
        {
          PyErr_SetString(PyExc_RuntimeError, "Error while computing symmetry operations.\n");
          bp::throw_error_already_set();
        }
        return result;
      }
    void expose_symmetry_operator()
    {
      bp::def
      (
        "get_space_group_symmetries",
        &get_space_group_symmetries,
        (
          bp::arg("lattice"),
          bp::arg("tolerance") = -1.0
        ),
        "Takes a (primitive) lattice and returns an array of symmetry "
        "operations for which the lattice is invariant.\n"
      );

      bp::def
      (
        "get_point_group_symmetries",
        &get_point_group_symmetries,
        (
          bp::arg("cell"),
          bp::arg("tolerance") = -1.0
        ),
        "Takes a cell and returns an array of symmetry "
        "operations for which the Bravais lattice is invariant.\n"
      );

      bp::scope scope = bp::class_<Crystal::SymmetryOperator>
        ( "SymmetryOperator", "SymmetryOperator" )
        .def(bp::init<atat::rMatrix3d const&>())
        .def(bp::init<atat::rVector3d const&>())
        .def(bp::init<atat::rMatrix3d const&, atat::rVector3d const&>())
        .def_readwrite("op", &Crystal::SymmetryOperator::op)
        .def_readwrite("trans", &Crystal::SymmetryOperator::trans)
        .def("invariant", &Crystal::SymmetryOperator::invariant, 
             (bp::arg("matrix"), bp::arg("tolerance")=types::tolerance),
             "Returns true if the matrix is invariant through this rotation.")
        .def("__call__", &Crystal::SymmetryOperator::operator())
        .def("__str__", &tostream<Crystal::SymmetryOperator>);

      expose_vector<Crystal::SymmetryOperator>
         ("Array", "An array of crystal.SymmetryOperator");
      bp::register_ptr_to_python< boost::shared_ptr< std::vector<Crystal::SymmetryOperator> > >();
    }

  }
} // namespace LaDa
