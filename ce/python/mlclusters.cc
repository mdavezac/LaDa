//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <iostream>
#include <sstream>
#include <complex>

#include <boost/python/class.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/default_call_policies.hpp>

#include <python/misc.hpp>
#include <python/std_vector.hpp>

#include "mlclusters.hpp"
#include "../mlclusters.h"


namespace LaDa
{
  namespace Python
  {
    CE::MLClusters::const_reference getvecitem( const CE::MLClusters& _vec, types::t_int _i )
    {
      const types::t_int dim(  _vec.size() );
      if( _i > dim or _i < -dim )
      {
        std::ostringstream sstr;
        throw std::out_of_range( "CE::MLClusters" );
      }
      return _vec[ size_t( _i < 0 ? dim + _i: _i ) ];
    }
    void setvecitem( CE::MLClusters& _vec, types::t_int _i,
                     CE::MLClusters::const_reference _a )
    {
      const types::t_int dim(  _vec.size() );
      if( _i > dim or _i < -dim )
      {
        std::ostringstream sstr;
        throw std::out_of_range( "atat vector" );
      }
      _vec[ size_t( _i < 0 ? types::t_int(dim) + _i: _i ) ] = _a;
    }

    boost::shared_ptr<CE::MLClusters> init( CE::MLCluster const &_cls )
    {
      boost::shared_ptr<CE::MLClusters> result(new CE::MLClusters);
      result->init(_cls);
      return result;
    }
    
    std::ostream& operator<<( std::ostream &_stream, CE::t_MLClusterClasses const &_cls )
    {
      CE::t_MLClusterClasses::const_iterator i_first = _cls.begin();
      CE::t_MLClusterClasses::const_iterator const i_end = _cls.end();
      for(; i_first != i_end; ++i_first) _stream << *i_first;
      return _stream;
    }

    boost::shared_ptr<CE::t_MLClusterClasses> init2( Crystal::Lattice const &_lat,
                                                     std::string const &_path,
                                                     std::string const &_genes )
      { return CE::read_clusters( _lat, _path, _genes); }


    void expose_mlclusters()
    {
      namespace bp = boost::python;
      expose_vector<CE::t_MLClusterClasses>
        ("MLClusterClasses", "An array of arrays of (presumably) equivalent multi-lattice clusters.\n")
        .def
        (
          "__init__",
          bp::make_constructor
          ( 
            &init2, 
            bp::default_call_policies(), 
            (bp::arg("lattice"), bp::arg("file"), bp::arg("genes")="")
          ),
          "Construct class of equivalent clusters from file. \n"
          "If genes is not empty, then it must be a string of 0 and 1, \n"
          "where 0 means that cluster in the input file will be ignored.\n"
        );
      
      bp::scope scope = bp::class_<CE::MLClusters>
        ("MLClusters", "An array of equivalent multi-lattice cluster (figure).")
          .def(bp::init<CE::MLClusters const&>())
          .def("__init__", bp::make_constructor(&init), "Creates a class of equivalent clusters from input.\n" )
          .add_property("eci", &CE::MLClusters::eci, "Interaction energy.")
          .def("order", &CE::MLClusters::order, "Returns the order of the cluster.\n")
          .def("__len__", &CE::MLClusters::size, "Returns the order of the cluster.\n")
          .def("__str__", &tostream<CE::MLClusters>, "Returns the order of the cluster.\n")
          .def("__getitem__", &getvecitem, bp::return_internal_reference<>())
          .def("__setitem__", &setvecitem);
      bp::register_ptr_to_python< boost::shared_ptr<CE::MLClusters> >();
      bp::register_ptr_to_python< boost::shared_ptr<CE::t_MLClusterClasses> >();
    }

  }
} // namespace LaDa
