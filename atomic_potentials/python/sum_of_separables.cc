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

    void expose_sumofseps()
    {
      namespace bp = boost::python;

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
        .def("append", &SumOfSeparables::push_back)
        .def("__iter__", &create_sumofsepsiter);

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

    }
  }
} // namespace LaDa
