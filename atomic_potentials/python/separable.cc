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
   
    struct SeparableIter
    {
      SeparableIter   ( atomic_potential::Separable &_sep )
                    : first_(true), cit_(_sep.begin()), cit_end_(_sep.end()) {}
      SeparableIter   ( SeparableIter const &_c )
                    : cit_(_c.cit_), cit_end_(_c.cit_end_), first_(_c.first_) {}

      SeparableIter &iter()  { return *this; }
      atomic_potential::Separable::iterator::reference next()
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
        return *cit_;
      }

      atomic_potential::Separable::iterator cit_;
      atomic_potential::Separable::iterator cit_end_;
      bool first_;
    };
    SeparableIter create_separableiter( atomic_potential::Separable & _ss )
      { return SeparableIter(_ss); }

    void expose_separable()
    {
      namespace bp = boost::python;

      bp::scope scope = bp::class_<Separable>("Separable", "A separables function.")
        .def(bp::init<Separable const&>())
        .def("__call__", &Separable::operator()<atomic_potential::VariableSet::t_Variables>)
        .def("__call__", &call_str<Separable>, bp::arg("structure"))
        .def("__call__", &call_rep<Separable>, bp::arg("representation"))
        .def("__len__", &Separable::size)
        .def("append", &Separable::push_back)
        .def("__iter__", &create_separableiter);
      
      bp::class_<SeparableIter>
      (
        "iterator", 
        "An iterator to separable functions.",
        bp::init<SeparableIter const&>()
      ).def("__iter__", &SeparableIter::iter, bp::return_internal_reference<1>() )
       .def("next", &SeparableIter::next, bp::return_internal_reference<1>() );
    }
  }
} // namespace LaDa
