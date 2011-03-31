#include "LaDaConfig.h"

#include <iostream>
#include <sstream>
#include <complex>

#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

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
      if( _i >= dim or _i <= -dim )
      {
        PyErr_SetString(PyExc_IndexError, "CE::MLCluster");
        boost::python::throw_error_already_set();
        static CE::MLCluster::value_type val;
        return val;
      }
      return _vec[ size_t( _i < 0 ? dim + _i: _i ) ];
    }
    void setvecitem( CE::MLCluster& _vec, types::t_int _i,
                     CE::MLCluster::const_reference _a )
    {
      const types::t_int dim(  _vec.size() );
      if( _i >= dim or _i <= -dim )
      {
        PyErr_SetString(PyExc_IndexError, "CE::MLCluster");
        boost::python::throw_error_already_set();
        return;
      }
      _vec[ size_t( _i < 0 ? types::t_int(dim) + _i: _i ) ] = _a;
    }

    void append_spin( CE::MLCluster& _self, bp::object _object)
    {
      bool found = false;
      CE::MLCluster::Spin spin;
      try { spin = bp::extract<CE::MLCluster::Spin>(_object); found = true;}
      catch(...) { PyErr_Clear(); };
      if(bp::len(_object) == 2)
        try
        { 
          spin.pos[0] = bp::extract<types::t_real>(_object[0][0]);
          spin.pos[1] = bp::extract<types::t_real>(_object[0][1]);
          spin.pos[2] = bp::extract<types::t_real>(_object[0][2]);
          spin.site = bp::extract<size_t>(_object[1]);
          found = true;
        }
        catch(...){ PyErr_Clear(); }
      if( (not found) and bp::len(_object) == 3 )
        try
        {
          spin.pos[0] = bp::extract<types::t_real>(_object[0]);
          spin.pos[1] = bp::extract<types::t_real>(_object[1]);
          spin.pos[2] = bp::extract<types::t_real>(_object[2]);
          spin.site = 0;
          found = true;
        }
        catch(...) { PyErr_Clear(); }
      if( (not found) and bp::len(_object) == 1 )
        try
        {
          spin.pos = bp::extract<math::rVector3d>(_object[0]);
          spin.site = 0;
          found = true;
        }
        catch(...) { PyErr_Clear(); }
      if( not found) 
      {
        PyErr_SetString(PyExc_ValueError, "Could not convert input to cluster spin.");
        bp::throw_error_already_set();
        return;
      }
      _self.push_back(spin);
    }
    void append_two( CE::MLCluster& _self, bp::list const &_pos, int _i )
    {
      if(bp::len(_pos) != 3)
      {
        PyErr_SetString(PyExc_ValueError, "Wrong dimension for vector.");
        bp::throw_error_already_set();
        return;
      }
      CE::MLCluster::Spin spin;
      spin.site = _i;
      spin.pos[0] = bp::extract<types::t_real>(_pos[0]);
      spin.pos[1] = bp::extract<types::t_real>(_pos[1]);
      spin.pos[2] = bp::extract<types::t_real>(_pos[2]);
      _self.push_back(spin);
    }

    CE::MLCluster::Spin pop(CE::MLCluster  &_self, types::t_int _i)
    {
      if(_i < 0) _i += _self.size();
      if(_i < 0 or _i >= _self.size())
      {
        PyErr_SetString(PyExc_IndexError, "Index out-of-range on CE cluster.");
        bp::throw_error_already_set();
        return CE::MLCluster::Spin();
      }
      CE::MLCluster :: iterator i_erase = _self.begin() + _i;
      CE::MLCluster :: Spin result = *i_erase;
      _self.erase(i_erase);
      return result;
    }

    void expose_mlcluster()
    {
      namespace bp = boost::python;
      bp::scope scope = bp::class_<CE::MLCluster>("MLCluster", "A multi-lattice cluster (figure).")
          .def(bp::init<CE::MLCluster const&>())
          .def_readwrite("origin", &CE::MLCluster::origin, "Cluster origin.")
          .def("apply_symmetry", &apply_symmetry, "Returns a transformed cluster.\n")
          .def("order", &CE::MLCluster::order, "Returns the order of the cluster.\n")
          .def("__str__", &tostream<CE::MLCluster>, "Returns the order of the cluster.\n")
          .def("__len__", &CE::MLCluster::size, "Returns the order of the cluster.\n")
          .def("__getitem__", &getvecitem, bp::return_internal_reference<>() )
          .def("__setitem__", &setvecitem)
          .def("pop", &pop)
          .def("clear", &CE::MLCluster::clear)
          .def("append", &append_spin)
          .def("append", &append_two,
               "Appends spin to current cluster.\n\n"
               "Input can consist of a `Spin`, a position and a site, "
               "or a only the position (site defaults to 0)." )
          .def_pickle( Python::pickle<CE::MLCluster>() );

      bp::class_<CE::MLCluster::Spin>("Spin", "A spin.")
        .def_readwrite("site", &CE::MLCluster::Spin::site, "Site index in lattice.")
        .add_property
        (
          "pos",
          make_getter(&CE::MLCluster::Spin::pos, bp::return_value_policy<bp::return_by_value>()),
          make_setter(&CE::MLCluster::Spin::pos, bp::return_value_policy<bp::return_by_value>()),
          "Relative position from origin. Numpy vector. Cartesian units." 
        );
    }

  }
} // namespace LaDa
