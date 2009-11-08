//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <iostream>
#include <sstream>
#include <complex>

#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/return_internal_reference.hpp>

#include "mlcluster.hpp"
#include "../mlcluster.h"

#include <python/misc.hpp>


namespace LaDa
{
  namespace Python
  {
    CE::MLCluster apply_symmetry( CE::MLCluster const &_self, Crystal::SymmetryOperator const &_op )
    {
      CE::MLCluster result(_self);
      result.apply_symmetry(_op);
      return result;
    }

    CE::MLCluster::const_reference getvecitem( const CE::MLCluster& _vec, types::t_int _i )
    {
      const types::t_int dim(  _vec.size() );
      if( _i > dim or _i < -dim )
      {
        std::ostringstream sstr;
        throw std::out_of_range( "CE::MLCluster" );
      }
      return _vec[ size_t( _i < 0 ? dim + _i: _i ) ];
    }
    void setvecitem( CE::MLCluster& _vec, types::t_int _i,
                     CE::MLCluster::const_reference _a )
    {
      const types::t_int dim(  _vec.size() );
      if( _i > dim or _i < -dim )
      {
        std::ostringstream sstr;
        throw std::out_of_range( "atat vector" );
      }
      _vec[ size_t( _i < 0 ? types::t_int(dim) + _i: _i ) ] = _a;
    }

    void expose_mlcluster()
    {
      namespace bp = boost::python;
      bp::scope scope = bp::class_<CE::MLCluster>("MLCluster", "A multi-lattice cluster (figure).")
          .def(bp::init<CE::MLCluster const&>())
          .add_property("origin", &CE::MLCluster::origin, "Cluster origin.")
          .def("apply_symmetry", &apply_symmetry, "Returns a transformed cluster.\n")
          .def("order", &CE::MLCluster::order, "Returns the order of the cluster.\n")
          .def("__str__", &tostream<CE::MLCluster>, "Returns the order of the cluster.\n")
          .def("__len__", &CE::MLCluster::size, "Returns the order of the cluster.\n")
          .def("__getitem__", &getvecitem, bp::return_internal_reference<>() )
          .def("__setitem__", &setvecitem);
      bp::class_<CE::MLCluster::Spin>("Spin", "A spin.")
        .add_property("site", &CE::MLCluster::Spin::site, "Site index in lattice.")
        .add_property("pos", &CE::MLCluster::Spin::pos, "Relative position from origin.");
    }

  }
} // namespace LaDa
