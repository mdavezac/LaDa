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
#include <boost/python/dict.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/with_custodian_and_ward.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <crystal/structure.h>
#include <python/misc.hpp>

#include "../representation.h"
#include "../functions.h"


namespace LaDa
{
  namespace Python
  {
    namespace bp = boost::python;
    typedef atomic_potential::Functions Functions;

    size_t getN () { return Functions::N; }

    struct Function
    {
      typedef Functions::arg_type::first_type arg_type;
      typedef Functions::result_type result_type;
      Function( bp::object const& _c ) : object(_c) {}
      Function( Function const& _c ) : object(_c.object) {}
      result_type operator()( arg_type const &_arg )
        { return bp::extract<result_type>( object(_arg) ); }
      bp::object object;
    };

    bool extract_tuple( bp::tuple const &_tuple,
                        Functions::t_Coefficient _coefs[] )
    {
      if( not bp::len( _tuple ) == Functions::N )
      {
        PyErr_SetString
        (
          PyExc_TypeError, 
          "Tuple is not of correct size.\n" 
        );
        return false;
      }
      try
      {
        for(size_t i(0); i < Functions::N; ++i)
          _coefs[i] = bp::extract<Functions::t_Coefficient>( _tuple[i] );
      }
      catch(...)
      {
        PyErr_SetString
        (
          PyExc_TypeError, 
          "Could not extract coefficients.\n"
        );
        return false;
      }
      return true;
    }    
    void  push_back_python_callable( Functions &_func, 
                                     bp::object const &_object, 
                                     bp::tuple const &_tuple )
    {    
      Functions::t_Coefficient coefs[Functions::N];
      if( not extract_tuple(_tuple, coefs)) { bp::throw_error_already_set(); return; }
      
      bp::dict const dict( _object.attr("__dict__") );
      if( not dict.has_key("__call__") )
      {
        PyErr_SetString
        (
          PyExc_TypeError, 
          "Object is not callable.\n" 
        );
        bp::throw_error_already_set();
        return;
      }
      _func.push_back( Function(_object), coefs);
    }

    template<class T_TYPE>
      void push_back_pow( Functions &_func, T_TYPE _pow, 
                          bp::tuple const &_tuple )
      {
        namespace bl = boost::lambda;

        Functions::t_Coefficient coefs[Functions::N];
        if( not extract_tuple(_tuple, coefs)) { bp::throw_error_already_set(); return; }

        typedef atomic_potential::numeric_type numeric_type;
        numeric_type (*ptr_func)(numeric_type, T_TYPE) = &std::pow;
        _func.push_back( bl::bind(ptr_func, bl::_1, bl::constant(_pow)), coefs);
      }

    atomic_potential::numeric_type constant( Functions::arg_type::first_type const& )
      { return 1e0; }

    void push_back_constant(Functions &_func, bp::tuple const &_tuple)
    {

      Functions::t_Coefficient coefs[Functions::N];
      if( not extract_tuple(_tuple, coefs)) { bp::throw_error_already_set(); return; }

      _func.push_back( &constant, coefs);
    }

    typedef boost::tuples::tuple< Functions::iterator, 
                                  Functions::iterator,
                                  bool > t_IterTuple;
    t_IterTuple& iter_self( t_IterTuple & _this ) { return _this; }
    t_IterTuple iter( Functions & _this )
      { return t_IterTuple(_this.begin(), _this.end(), true); }
    Functions::iterator::reference next( t_IterTuple & _this )
    {
      namespace bt = boost::tuples;
      if( bt::get<2>(_this) ) bt::get<2>(_this) = false;
      else if( bt::get<0>(_this) != bt::get<1>(_this) ) ++bt::get<0>(_this);
      if( bt::get<0>(_this) == bt::get<1>(_this) )
      {
        PyErr_SetString( PyExc_StopIteration, "End-of-range.\n");
        bp::throw_error_already_set();
        return Functions::iterator::reference();
      }
      return *bt::get<0>(_this);
    }

    Functions::result_type call_1( Functions::iterator::reference const &_ref,
                                   Functions::arg_type::first_type _a )
      { return _ref(_a); }
    Functions::result_type call_2( Functions::iterator::reference const &_ref,
                                   Functions::arg_type const &_a )
      { return _ref(_a); }
    Functions::result_type call_3( Functions::iterator::reference const &_ref,
                                   bp::tuple const &_a )
    {
      Functions::arg_type a;
      if( not bp::len( _a ) == 2 )
      {
        PyErr_SetString(PyExc_TypeError, "Tuple is not a 2-tuple.\n");
        bp::throw_error_already_set();
        return -1;
      }
      try
      {
        return _ref
               ( 
                 Functions::arg_type
                 (
                   bp::extract<Functions::arg_type::first_type>( _a[0] ),
                   bp::extract<Functions::arg_type::second_type>( _a[1] )
                 )
               );
      }
      catch( ... )
      {
        PyErr_SetString
        (
          PyExc_TypeError, 
          "Could not extract variable.\n"
        );
        return -1;
      }
    }


    atomic_potential::numeric_type getitem( Functions::iterator::reference const &_ref,
                                            size_t _i )
      {
        if( _i >= Functions::N )
        {
          PyErr_SetString( PyExc_IndexError, "Coefficient index out of range.\n");
          bp::throw_error_already_set();
          static atomic_potential::numeric_type a(-1);
          return a;
        }
        return _ref[_i];
      }
    atomic_potential::numeric_type setitem( Functions::iterator::reference &_ref,
                                            size_t _i, atomic_potential::numeric_type _a )
      {
        if( _i >= Functions::N )
        {
          PyErr_SetString( PyExc_IndexError, "Coefficient index out of range.\n");
          bp::throw_error_already_set();
          static atomic_potential::numeric_type a(-1);
          return a;
        }
        _ref[_i] = _a;
      }

    size_t len(Functions::iterator::reference const& ) { return Functions::N; }
    
    void expose_functions()
    {
      bp::scope scope = bp::class_<Functions>("Functions", "A function of a single coordinate.")
        .def(bp::init<Functions const&>())
        .def("__call__", &Functions::operator())
        .def("__len__", &Functions::size)
        .def("append", &push_back_python_callable )
        .def("append_pow", &push_back_pow<atomic_potential::numeric_type>)
        .def("append_pow", &push_back_pow<int>)
        .def("append_constant", &push_back_constant)
        .def("__str__", &tostream<Functions>)
        .def("__iter__", &iter, bp::with_custodian_and_ward_postcall<1,0>())
        .add_static_property("N", &getN);

      bp::class_<t_IterTuple>
      (
        "iterator", 
        "An iterator to functions.",
        bp::no_init
      ).def("__iter__", &iter_self, bp::return_internal_reference<1>())
       .def("next", &next, bp::with_custodian_and_ward_postcall<1,0>());
 
      bp::class_<Functions::iterator::reference>
      (
        "reference", 
        "reference to a function obtained from an iterator.",
        bp::init<Functions::iterator::reference const&>()
      ).def("__getitem__", &getitem)
       .def("__setitem__", &setitem)
       .def("__len__", &len)
       .def("__call__", &call_1, "Calls function with a single real argument. Coefs are ignored.")
       .def("__call__", &call_2, "Calls function with a real and a variable argument.")
       .def("__call__", &call_3, "Calls function with a tuple of a real and a variable argument.");
 
    }
  }
} // namespace LaDa
