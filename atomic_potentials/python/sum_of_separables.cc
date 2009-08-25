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

#include <crystal/structure.h>

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

    atomic_potential::numeric_type constant( Functions::arg_type::first_type const& ) { return 1e0; }

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
        .def("append", &push_back_python_callable )
        .def("append_pow", &push_back_pow<atomic_potential::numeric_type>)
        .def("append_pow", &push_back_pow<int>)
        .def("append_constant", &push_back_constant);
      
      bp::class_<Separable>("Separable", "A separables function.")
        .def(bp::init<Separable const&>())
        .def("__call__", &Separable::operator()<atomic_potential::VariableSet::t_Variables>)
        .def("__call__", &call_str<Separable>, bp::arg("structure"))
        .def("__call__", &call_rep<Separable>, bp::arg("representation"))
        .def("append", &Separable::push_back);
      
      bp::class_<SumOfSeparables>("SumOfSeparables", "A sum of separables function.")
        .def(bp::init<SumOfSeparables const&>())
        .def("__call__", &SumOfSeparables::operator()<atomic_potential::VariableSet::t_Variables>)
        .def("__call__", &call_str<SumOfSeparables>, bp::arg("structure"))
        .def("__call__", &call_rep<SumOfSeparables>, bp::arg("representation"))
        .def("append", &SumOfSeparables::push_back);
    }
  }
} // namespace LaDa
