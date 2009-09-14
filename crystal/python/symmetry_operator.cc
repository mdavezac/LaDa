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
    void expose_symmetry_operator()
    {
      bp::def
      (
        "get_symmetries",
        &Crystal::get_symmetries,
        "Takes a lattice and returns an array of symmetry operations for which the lattice is invariant.\n"
      );

      bp::scope scope = bp::class_<Crystal::SymmetryOperator>( "SymmetryOperator", "SymmetryOperator" )
        .def(bp::init<atat::rMatrix3d const&>())
        .def(bp::init<atat::rVector3d const&>())
        .def(bp::init<atat::rMatrix3d const&, atat::rVector3d const&>())
        .def_readwrite("op", &Crystal::SymmetryOperator::op)
        .def_readwrite("trans", &Crystal::SymmetryOperator::trans)
        .def("__call__", &Crystal::SymmetryOperator::operator())
        .def("__str__", &tostream<Crystal::SymmetryOperator>);

      expose_vector<Crystal::SymmetryOperator>
         ("Array", "An array of crystal.SymmetryOperator");
      bp::register_ptr_to_python< boost::shared_ptr< std::vector<Crystal::SymmetryOperator> > >();
    }

  }
} // namespace LaDa
