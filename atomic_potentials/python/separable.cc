//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <cmath>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/object.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <crystal/structure.h>
#include <python/misc.hpp>

#include "../representation.h"
#include "../separable.h"


namespace LaDa
{
  namespace Python
  {
    namespace bp = boost::python;
    typedef atomic_potential::Separable Separable;

    template<class T_TYPE> 
      typename T_TYPE::result_type call_rep( T_TYPE const& _sep,
                                             atomic_potential::Representation const & _rep )
      {
        typename T_TYPE::result_type result(0);
        atomic_potential::Representation::const_iterator i_first = _rep.begin();
        atomic_potential::Representation::const_iterator const i_end = _rep.end();
      
        for(; i_first != i_end; ++i_first)
          result += _sep(i_first->variables) * i_first->weight;
        return result;
      }
    template<class T_TYPE>
      typename T_TYPE::result_type call_str( T_TYPE const& _sep,
                                             Crystal::TStructure<std::string> const & _str )
      {
        atomic_potential::Representation const representation(_str, _sep.nb_coordinates() );
        return call_rep(_sep, representation);
      }

    typedef boost::tuples::tuple< Separable::iterator, 
                                  Separable::iterator,
                                  bool > t_IterTuple;
    t_IterTuple& iter_self( t_IterTuple & _this ) { return _this; }
    t_IterTuple iter( Separable & _this )
      { return t_IterTuple(_this.begin(), _this.end(), true); }
    Separable::iterator::reference next( t_IterTuple & _this )
    {
      namespace bt = boost::tuples;
      if( bt::get<2>(_this) ) bt::get<2>(_this) = false;
      else if( bt::get<0>(_this) != bt::get<1>(_this) ) ++bt::get<0>(_this);
      if( bt::get<0>(_this) == bt::get<1>(_this) )
      {
        PyErr_SetString( PyExc_StopIteration, "End-of-range.\n");
        bp::throw_error_already_set();
        --bt::get<0>(_this);
      }
      return *bt::get<0>(_this);
    }
   
    void expose_separable()
    {
      bp::scope scope = bp::class_<Separable>("Separable", "A separables function.")
        .def(bp::init<Separable const&>())
        .def("__call__", &Separable::operator()<atomic_potential::VariableSet::t_Variables>)
        .def("__call__", &call_str<Separable>, bp::arg("structure"))
        .def("__call__", &call_rep<Separable>, bp::arg("representation"))
        .def("__len__", &Separable::size)
        .def("append", &Separable::push_back)
        .def("__iter__", &iter, bp::with_custodian_and_ward_postcall<1,0>());
      
      bp::class_<t_IterTuple>
      (
        "iterator", 
        "An iterator to separable functions.",
        bp::no_init
      ).def("__iter__", &iter_self, bp::return_internal_reference<1>() )
       .def("next", &next, bp::return_internal_reference<1>() );
    }
  }
} // namespace LaDa
