//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <cmath>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/tuple/tuple.hpp>

#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/object.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <crystal/structure.h>
#include <python/misc.hpp>

#include "../representation.h"
#include "../sum_of_separables.h"


namespace LaDa
{
  namespace Python
  {
    namespace bp = boost::python;

    typedef LaDa::atomic_potential::SumOfSeparables SumOfSeparables;
    typedef SumOfSeparables::t_Function Separable;
    typedef Separable::t_Function Functions;
    typedef atomic_potential::Representation Representation;
    typedef atomic_potential::VariableSet VariableSet;

    template<class T_TYPE> 
      typename T_TYPE::result_type call_rep( T_TYPE const& _sep,
                                             Representation const & _rep )
      {
        typename T_TYPE::result_type result(0);
        Representation::const_iterator i_first = _rep.begin();
        Representation::const_iterator const i_end = _rep.end();
      
        for(; i_first != i_end; ++i_first)
          result += _sep(i_first->variables) * i_first->weight;
        return result;
      }
    template<class T_TYPE>
      typename T_TYPE::result_type call_str( T_TYPE const& _sep,
                                             Crystal::TStructure<std::string> const & _str )
      {
        Representation const representation(_str, (_sep.nb_coordinates()+2)/3);
        types::t_real const result( call_rep(_sep, representation) );
        return result;
      }
   
    typedef boost::tuples::tuple< SumOfSeparables::iterator, 
                                  SumOfSeparables::iterator,
                                  bool > t_IterTuple;
    t_IterTuple& iter_self( t_IterTuple & _this ) { return _this; }
    t_IterTuple iter( SumOfSeparables & _this )
      { return t_IterTuple(_this.begin(), _this.end(), true); }
    boost::python::tuple next( t_IterTuple & _this )
    {
      namespace bt = boost::tuples;
      if( bt::get<2>(_this) ) bt::get<2>(_this) = false;
      else if( bt::get<0>(_this) != bt::get<1>(_this) ) ++bt::get<0>(_this);
      if( bt::get<0>(_this) == bt::get<1>(_this) )
      {
        PyErr_SetString( PyExc_StopIteration, "End-of-range.\n");
        bp::throw_error_already_set();
        return bp::make_tuple(-1);
      }
      return bp::make_tuple( bt::get<0>(*bt::get<0>(_this)), 
                             bt::get<1>(*bt::get<0>(_this)) );
    }

    void expose_sumofseps()
    {
      bp::scope scope = bp::class_<SumOfSeparables>
              ("SumOfSeparables", "A sum of separables function.")
        .def(bp::init<SumOfSeparables const&>())
        .def("__call__", &SumOfSeparables::operator()<atomic_potential::VariableSet::t_Variables>)
        .def("__call__", &call_str<SumOfSeparables>, bp::arg("structure"))
        .def("__call__", &call_rep<SumOfSeparables>, bp::arg("representation"))
        .add_property
        (
          "nb_coordinates", 
          &SumOfSeparables::nb_coordinates,
          "Number of coordinates.\n" 
        )
        .def("__len__", &SumOfSeparables::size)
        .def("append", &SumOfSeparables::push_back)
        .def("__iter__", &iter, bp::with_custodian_and_ward_postcall<1,0>());

      bp::class_<t_IterTuple>
      (
        "SumOfSeparablesIterator", 
        "An iterator to separable functions.",
        bp::no_init
      ).def("__iter__", &iter_self, bp::return_internal_reference<1>() )
       .def("next", &next, bp::with_custodian_and_ward_postcall<1,0>());

    }
  }
} // namespace LaDa
