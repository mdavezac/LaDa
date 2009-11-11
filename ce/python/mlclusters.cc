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
#include <boost/python/list.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/default_call_policies.hpp>

#include <pyublas/numpy.hpp>

#include <python/misc.hpp>

#include "mlclusters.hpp"
#include "../mlclusters.h"
#include "../find_pis.h"


namespace LaDa
{
  namespace Python
  {
    CE::MLClusters::const_reference getvecitem( const CE::MLClusters& _vec, types::t_int _i )
    {
      const types::t_int dim(  _vec.size() );
      if( _i >= dim or _i <= -dim )
      {
        PyErr_SetString(PyExc_IndexError, "CE::MLCluster");
        boost::python::throw_error_already_set();
        static CE::MLClusters::value_type val;
        return val;
      }
      return _vec[ size_t( _i < 0 ? dim + _i: _i ) ];
    }
    void setvecitem( CE::MLClusters& _vec, types::t_int _i,
                     CE::MLClusters::const_reference _a )
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
    CE::t_MLClusterClasses::const_reference getvecitem2( const CE::t_MLClusterClasses& _vec, types::t_int _i )
    {
      const types::t_int dim(  _vec.size() );
      if( _i >= dim or _i <= -dim )
      {
        PyErr_SetString(PyExc_IndexError, "CE::MLCluster");
        boost::python::throw_error_already_set();
        static CE::t_MLClusterClasses::value_type val;
        return val;
      }
      return _vec[ size_t( _i < 0 ? dim + _i: _i ) ];
    }
    void setvecitem2( CE::t_MLClusterClasses& _vec, types::t_int _i,
                      CE::t_MLClusterClasses::const_reference _a )
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

    types::t_real call( CE::t_MLClusterClasses const &_cls, Crystal::Structure const &_str )
    {
      std::vector<types::t_real> pis;
      CE::find_pis(_cls, _str, pis);
      types::t_real result(0);
      LADA_ASSERT( _cls.size() == pis.size(), "Inconsistent sizes.\n")
      std::vector<types::t_real> :: const_iterator i_pi = pis.begin();
      std::vector<types::t_real> :: const_iterator const i_pi_end = pis.end();
      CE::t_MLClusterClasses::const_iterator i_class = _cls.begin();
      for(; i_pi != i_pi_end; ++i_pi, ++i_class)
        result += (*i_pi) * i_class->eci;
      return result;
    }
     
    pyublas::numpy_vector<types::t_real>  pis( CE::t_MLClusterClasses const &_cls,
                                                  Crystal::Structure const &_str )
    {
      pyublas::numpy_vector<types::t_real> pis_;
      CE::find_pis(_cls, _str, pis_);
      return pis_;
    }

    boost::shared_ptr<CE::t_MLClusterClasses> init2( Crystal::Lattice const &_lat,
                                                     std::string const &_path,
                                                     std::string const &_genes )
      { return CE::read_clusters( _lat, _path, _genes); }
    boost::shared_ptr<CE::t_MLClusterClasses> init3(std::string const &_path, bool _b)
      { return CE::load_mlclusters(_path, _b); }

    boost::shared_ptr<CE::t_MLClusterClasses> copy_constructor( CE::t_MLClusterClasses const & _ob )
      { return boost::shared_ptr<CE::t_MLClusterClasses>( new CE::t_MLClusterClasses(_ob) ); }
    boost::shared_ptr<CE::t_MLClusterClasses> object_constructor( const boost::python::list& _ob )
    {
      namespace bp = boost::python;
      const size_t n( bp::len( _ob ) );
      boost::shared_ptr<CE::t_MLClusterClasses> result( new CE::t_MLClusterClasses );
      result->reserve( n );
      for( size_t i(0); i < n; ++i )
        result->push_back( bp::extract<CE::MLClusters>( _ob[i] ) );
      return result;
    }

    void expose_mlclusters()
    {
      namespace bp = boost::python;
      bp::register_ptr_to_python< boost::shared_ptr<CE::t_MLClusterClasses> >();
      bp::class_<CE::t_MLClusterClasses>
        ("MLClusterClasses", "An array of arrays of (presumably) equivalent multi-lattice clusters.\n")
        .def( "__init__", bp::make_constructor( &copy_constructor ) )
        .def( "__init__", bp::make_constructor( &object_constructor ) )
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
        )
        .def
        (
          "__init__",
          bp::make_constructor
          ( 
            &init3, 
            bp::default_call_policies(), 
            (bp::arg("file"), bp::arg("is_multi")=true)
          ),
          "Construct class of equivalent clusters from an xml file. \n"
          "If is_multi is false, then expects a single lattice format.\n" 
        )
        .def("__getitem__", &getvecitem2, bp::return_internal_reference<>())
        .def("__setitem__", &setvecitem2)
        .def("__str__", &tostream<CE::t_MLClusterClasses> )
        .def("clear", &CE::t_MLClusterClasses :: clear )
        .def("__len__", &CE::t_MLClusterClasses::size)
        .def("__call__", &call, "Returns energy of structure.\n")
        .def("pis", &pis, "Return numpy vector corresponding to structure pis.\n");
      
      bp::scope scope = bp::class_<CE::MLClusters>
        ("MLClusters", "An array of equivalent multi-lattice cluster (figure).")
          .def(bp::init<CE::MLClusters const&>())
          .def( "__init__", bp::make_constructor(&init),
                "Creates a class of equivalent clusters from input.\n" )
          .add_property("eci", &CE::MLClusters::eci, "Interaction energy.")
          .def("order", &CE::MLClusters::order, "Returns the order of the cluster.\n")
          .def("__len__", &CE::MLClusters::size, "Returns the number of equivalent cluster.\n")
          .def("__str__", &tostream<CE::MLClusters>, "Returns the order of the cluster.\n")
          .def("__getitem__", &getvecitem, bp::return_internal_reference<>())
          .def("__setitem__", &setvecitem);
      bp::register_ptr_to_python< boost::shared_ptr<CE::MLClusters> >();
    }

  }
} // namespace LaDa
