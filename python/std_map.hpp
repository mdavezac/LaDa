#ifndef _LADA_OPT_PYTHON_STD_MAP_HPP_
#define _LADA_OPT_PYTHON_STD_MAP_HPP_

#include "LaDaConfig.h"

#include <vector>
#include <string>
#include <sstream>
#include <boost/python/class.hpp>
#include <boost/python/list.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/manage_new_object.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/with_custodian_and_ward.hpp>

#include <opt/debug.h>

namespace LaDa
{
  namespace Python
  {
    namespace details
    {
#     include "std_map.impl.hpp"
    }

    template< class T_MAP, template< class > class T_DEREF >
      void expose_iterator( const std::string &_name )
      {
        typedef T_MAP t_Map;
        typedef details::map_iterator< t_Map, T_DEREF >  t_iterator;
        namespace bp = boost::python;
        bp::class_< t_iterator >( _name.c_str() )
          .def( bp::init< t_iterator >() )
          .def( "__iter__", &t_iterator::iter, bp::return_internal_reference<1>() )
          .def( "next", &t_iterator::next, bp::return_internal_reference<1>() )
          .def( "check", &t_iterator::check );
      }
    template< class T_MAP >
      void expose_map( const std::string &_name, const std::string &_docstring )
      {
        typedef T_MAP t_Map;
        typedef details::map_iterator< t_Map, details::deref_item >  t_itemiterator;
        typedef details::map_iterator< t_Map, details::deref_str_key >  t_keyiterator;
        namespace bp = boost::python;
        expose_iterator< t_Map, details::deref_item >( ("__itemiterator_" + _name).c_str() );
        expose_iterator< t_Map, details::deref_str_key >( ("__keyiterator_" + _name).c_str() );

        bp::class_< t_Map >( _name.c_str(), _docstring.c_str() )
          .def( "__init__",
                bp::make_constructor( details::default_constructor< T_MAP > ) )
          .def( "__init__",
                bp::make_constructor( details::copy_constructor< T_MAP > ) )
          .def( "__init__",
                bp::make_constructor( details::dict_constructor< T_MAP > ) )
          .def( "clear", &t_Map::clear )
          .def( "copy", &details::shallow_copy<T_MAP>, bp::return_internal_reference<1>() )
          .def( "__getitem__", &details::getitem<T_MAP> )
          .def( "__setitem__", &details::setitem<T_MAP> )
          .def( "get", &details::get<T_MAP> ) 
          .def( "get", &details::get_default<T_MAP> )
          .def( "setdefault", &details::set_default<T_MAP> )
          .def( "setdefault", &details::set_default_none<T_MAP> )
          .def( "__len__", &T_MAP::size )
          .def( "__contains__", &details::has_key<T_MAP> )
          .def( "has_key", &details::has_key<T_MAP> )
          .def( "items", &details::items<T_MAP>, bp::return_internal_reference<1>() )
          .def( "values", &details::values<T_MAP>, bp::return_internal_reference<1>() )
          .def( "keys", &details::keys<T_MAP>, bp::with_custodian_and_ward_postcall<1,0>() )
          .def( "pop", &details::pop<T_MAP> )
          .def( "popitem", &details::popitem<T_MAP> )
          .def( "__delitem__", &details::pop<T_MAP> )
          .def( "__iter__", &details::iter<t_itemiterator>, bp::return_internal_reference<1>() )
          .def( "iterkeys", &details::iter<t_keyiterator>, bp::return_internal_reference<1>() )
          .def( "__str__", &details::print<T_MAP> );
      }

    namespace details
    {
#     include "std_map.impl.hpp"
    }
  }
} // namespace LaDa
#endif
