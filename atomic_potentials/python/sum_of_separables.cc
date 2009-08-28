//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <cmath>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

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
    typedef LaDa::atomic_potential::SumOfSeparables SumOfSeparables;
    typedef SumOfSeparables::t_Function Separable;
    typedef Separable::t_Function Functions;

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
   
    struct Function
    {
      typedef Functions::arg_type::first_type arg_type;
      typedef Functions::result_type result_type;
      Function( boost::python::object const& _c ) : object(_c) {}
      Function( Function const& _c ) : object(_c.object) {}
      result_type operator()( arg_type const &_arg )
        { return boost::python::extract<result_type>( object(_arg) ); }
      boost::python::object object;
    };

    bool extract_tuple( boost::python::tuple const &_tuple,
                        Functions::t_Coefficient _coefs[Functions::N] )
    {
      namespace bp = boost::python;
      if( not bp::len( _tuple ) == Functions::N )
      {
        PyErr_SetString
        (
          PyExc_TypeError, 
          "Tuple is not a 2-tuple.\n" 
        );
        return false;
      }
      try
      {
        Functions::t_Coefficient coefs[Functions::N];
        for(size_t i(0); i < Functions::N; ++i)
          coefs[i] = bp::extract<Functions::t_Coefficient>( _tuple[i] );
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
                                     boost::python::object const &_object, 
                                     boost::python::tuple const &_tuple )
    {    
      namespace bp = boost::python;
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

    struct reference_ss
    {
      reference_ss   ( SumOfSeparables::iterator const &_it )
                   : sep_(_it->get<0>()), coef_(_it->get<1>()) {}
      reference_ss   ( reference_ss const & _ref )
                   : sep_(_ref.sep_), coef_(_ref.coef_) {}
      void set_sep( atomic_potential::Separable & _s ) { sep_ = _s; }
      atomic_potential::Separable const & get_sep() const { return sep_; }
      atomic_potential::numeric_type get_coef() const { return coef_; }
      void set_coef(atomic_potential::numeric_type _t) { coef_ = _t; }
      atomic_potential::Separable & sep_;
      atomic_potential::numeric_type & coef_;
    };

    struct SumOfSepsIter
    {
      SumOfSepsIter   ( atomic_potential::SumOfSeparables &_sumofseps )
                    : first_(true), cit_(_sumofseps.begin()), cit_end_(_sumofseps.end()) {}
      SumOfSepsIter   ( SumOfSepsIter const &_c )
                    : cit_(_c.cit_), cit_end_(_c.cit_end_), first_(_c.first_) {}

      SumOfSepsIter &iter()  { return *this; }
      reference_ss next()
      {
        namespace bp = boost::python;
        if( first_ ) first_ = false; 
        else 
        {
          ++cit_;
          if( cit_ == cit_end_ )
          {
            PyErr_SetString
            (
              PyExc_StopIteration, 
              "Error while computing transform to smith normal form.\n" 
            );
            bp::throw_error_already_set();
            --cit_;
          }
        }
        return cit_;
      }

      atomic_potential::SumOfSeparables::iterator cit_;
      atomic_potential::SumOfSeparables::iterator cit_end_;
      bool first_;
    };
    SumOfSepsIter create_sumofsepsiter( atomic_potential::SumOfSeparables & _ss )
      { return SumOfSepsIter(_ss); }

    template<class T_TYPE>
      void push_back_pow( Functions &_func, T_TYPE _pow, 
                          boost::python::tuple const &_tuple )
      {
        namespace bp = boost::python;
        namespace bl = boost::lambda;

        Functions::t_Coefficient coefs[Functions::N];
        if( not extract_tuple(_tuple, coefs)) { bp::throw_error_already_set(); return; }

        typedef atomic_potential::numeric_type numeric_type;
        numeric_type (*ptr_func)(numeric_type, T_TYPE) = &std::pow;
        _func.push_back( bl::bind(ptr_func, bl::_1, bl::constant(_pow)), coefs);
      }

    atomic_potential::numeric_type constant( Functions::arg_type::first_type const& )
      { return 1e0; }

    void push_back_constant(Functions &_func, boost::python::tuple const &_tuple)
    {
      namespace bp = boost::python;

      Functions::t_Coefficient coefs[Functions::N];
      if( not extract_tuple(_tuple, coefs)) { bp::throw_error_already_set(); return; }

      _func.push_back( &constant, coefs);
    }
    
    void expose_sumofseps()
    {
      namespace bp = boost::python;

      bp::class_<Functions>("Functions", "A function of a single coordinate.")
        .def(bp::init<Functions const&>())
        .def("__call__", &Functions::operator())
        .def("__len__", &Functions::size)
        .def("append", &push_back_python_callable )
        .def("append_pow", &push_back_pow<atomic_potential::numeric_type>)
        .def("append_pow", &push_back_pow<int>)
        .def("append_constant", &push_back_constant)
        .def("__str__", &tostream<Functions>);
      
      bp::class_<Separable>("Separable", "A separables function.")
        .def(bp::init<Separable const&>())
        .def("__call__", &Separable::operator()<atomic_potential::VariableSet::t_Variables>)
        .def("__call__", &call_str<Separable>, bp::arg("structure"))
        .def("__call__", &call_rep<Separable>, bp::arg("representation"))
        .def("__len__", &Separable::size)
        .def("append", &Separable::push_back);
      
      bp::class_<SumOfSepsIter>
      (
        "SumOfSeparablesIter", 
        "An iterator to separable functions.",
        bp::init<SumOfSepsIter const&>()
      ).def("__iter__", &SumOfSepsIter::iter, bp::return_internal_reference<1>() )
       .def("next", &SumOfSepsIter::next);

      bp::class_<reference_ss>
      (
        "SumOfSeparablesIterRef", 
        "Dereferenced iterator to separable functions.",
        bp::init<reference_ss const&>()
      ).add_property
       (
         "separables", 
         bp::make_function
         (
           &reference_ss::get_sep, 
           bp::return_internal_reference<>()
         ), 
         bp::make_function
         (
           &reference_ss::set_sep, 
           bp::return_internal_reference<>()
         ), 0
       ) 
       .add_property("coefficient", &reference_ss::get_coef, &reference_ss::set_coef);

      bp::class_<SumOfSeparables>("SumOfSeparables", "A sum of separables function.")
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
        .def("append", &SumOfSeparables::push_back)
        .def("__iter__", &create_sumofsepsiter);
    }
  }
} // namespace LaDa
