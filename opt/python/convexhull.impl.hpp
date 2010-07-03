#ifndef _PYTHONLADA_CONVEXHUL_IMPL_H_
#define _PYTHONLADA_CONVEXHUL_IMPL_H_

#include "LaDaConfig.h"

#include <boost/python/object.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/str.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/with_custodian_and_ward.hpp>
#include <boost/tuple/tuple.hpp>
#include <string>

#include "../convex_hull.h"
#include "../types.h"
#include "../debug.h"

namespace LaDa
{
  namespace Python
  {
    template< class T_TYPE > 
      struct PythonConvexHull
      {
        typedef LaDa::opt::ConvexHull::FakeObject< T_TYPE >   t_FakeObject;
        typedef LaDa::opt::ConvexHull::Base< t_FakeObject > t_CHBase;
        
        static void add( t_CHBase &_ch, T_TYPE const &_t, types::t_real _x,
                        types::t_real _energy )
        {
          t_FakeObject fake( _t );
          LaDa::opt::ConvexHull::NullObject::x = _x;
          _ch.add( _energy, fake );
        }
        
        typedef typename t_CHBase::const_iterator t_cit;
        typedef boost::tuples::tuple<t_cit, t_cit, bool> t_range;
        static t_range iter(t_CHBase const &_ch )
            { return t_range(_ch.begin_vertex(), _ch.end_vertex(), true); }
        static t_range& iter2(t_range &_this) { return _this; }
        static boost::python::tuple next(t_range &_this)
        {
          namespace bp = boost::python;
          namespace bt = boost::tuples;
          if( bt::get<2>(_this) ) bt::get<2>(_this) = false; 
          else 
          {
            ++bt::get<0>(_this);
            if( bt::get<0>(_this) == bt::get<1>(_this) )
            {
              PyErr_SetString(  PyExc_StopIteration, "End-of-range in convex-hull iteration." );
              bp::throw_error_already_set();
              --bt::get<0>(_this);
            }
          }
          return bp::make_tuple( bt::get<0>(_this)->object.object,
                                 bt::get<0>(_this)->x, 
                                 bt::get<0>(_this)->y );
        }
      };

    template< class T_TYPE >
      typename PythonConvexHull<T_TYPE>::t_CHBase* empty()
        { return new typename PythonConvexHull<T_TYPE>::t_CHBase(); }
    template< class T_TYPE >
      typename PythonConvexHull<T_TYPE>::t_CHBase* 
        copy( const typename PythonConvexHull<T_TYPE>::t_CHBase &_o )
          { return new typename PythonConvexHull<T_TYPE>::t_CHBase( _o ); }
   
  }

  namespace opt 
  { 

    namespace ConvexHull
    {
      template<>
        std::ostream& operator<<( std::ostream &_str,
                                  const FakeObject<boost::python::object> & _o )
        {
          boost::python::str str( _o.object );
          const std::string stringobject = boost::python::extract< const std::string>(str);
          return _str << std::setw(30) << stringobject;
        }
    }
   
  }

  namespace Python
  {
    template< class T_TYPE > void exposeConvexHull( const std::string &_name )
    {
      namespace bp = boost::python;
      bp::scope scope = 
        bp::class_< typename PythonConvexHull<T_TYPE>::t_CHBase >( _name.c_str() )
          .def( "__init__", bp::make_constructor( copy<T_TYPE> ) )
          .def( "__init__", bp::make_constructor( empty<T_TYPE> ) )
          .def( "__call__", &PythonConvexHull<T_TYPE>::t_CHBase::evaluate,
                "Evaluates convex-hull at C{x}." )
          .def( "add", &PythonConvexHull<T_TYPE>::add, 
                (bp::arg("item"), bp::arg("x"), bp::arg("y")),
                "Tries to add item to convex hull\n\n"
                "@param item: The item which is to be considered for a breaking point.\n"
                "@param x: concentration (x-axis) of the item.\n"
                "@param y: Energy (y-axis) of the item.\n" )
          .def("__iter__", &PythonConvexHull<T_TYPE>::iter, "Iterates over breaking-points.")
          .def( "__str__", &PythonConvexHull<T_TYPE>::t_CHBase::print );

      bp::class_< typename PythonConvexHull<T_TYPE>::t_range >( (_name+"_iter").c_str()  )
          .def("__iter__", &PythonConvexHull<T_TYPE>::iter2,
               bp::return_internal_reference<1>() )
          .def("next", &PythonConvexHull<T_TYPE>::next );
    }

  } //en of Python namespace
} // namespaecd LaDa

#endif
