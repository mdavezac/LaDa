#ifndef _LADA_OPT_PYTHON_STD_VECTOR_HPP_
#define _LADA_OPT_PYTHON_STD_VECTOR_HPP_

#include "LaDaConfig.h"

#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <boost/python/iterator.hpp>
#include <boost/python/class.hpp>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/with_custodian_and_ward.hpp>
#include <boost/python/copy_non_const_reference.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/return_arg.hpp>
#include <boost/python/str.hpp>

#include <root_exceptions.h>
#include <opt/debug.h>

#include "exceptions.h"

namespace LaDa
{
  namespace python
  {
    namespace bp = boost::python;
    namespace details
    {
      template<class T> std::auto_ptr< std::vector<T> > default_constructor()
        { return std::auto_ptr< std::vector<T> >(new std::vector<T>); }
      template<class T> std::auto_ptr< std::vector<T> > copy_constructor(std::vector<T> const &_c)
        { return std::auto_ptr< std::vector<T> >(new std::vector<T>(_c)); }
      template<class T>
        std::auto_ptr< std::vector<T> > constructor(bp::object _b)
        {
          std::auto_ptr< std::vector<T> > result(new std::vector<T>);
          if(not PyObject_HasAttrString(_b.ptr(), "__iter__"))
            PyException<error::AttributeError>::throw_error("Object does not implement __iter__.");
          bp::stl_input_iterator<T> i_first(_b);
          bp::stl_input_iterator<T> const i_end;
          for(; i_first != i_end; ++i_first) result->push_back(*i_first);
          return result;
        }
      template<class T>
        std::string __str__( const std::vector<T>& _ob )
        {
          if(_ob.size() == 0) return "[]";
          std::ostringstream sstr;
          typename std::vector<T>::const_iterator i_first = _ob.begin();
          typename std::vector<T>::const_iterator const i_end = _ob.end();
          sstr << "[ " << *i_first;
          for(++i_first; i_first != i_end; ++i_first)
            sstr << ", " << *i_first;
          sstr << "]";
          return sstr.str();
        }
      template<class T>
        void extend( std::vector<T> & _self, bp::object _o)
        {
          if(not PyObject_HasAttrString(_o.ptr(), "__iter__"))
            python::PyException<error::AttributeError>::throw_error("Object does not implement __iter__.");
          bp::stl_input_iterator<T> i_first(_o);
          bp::stl_input_iterator<T> const i_end;
          for(; i_first != i_end; ++i_first) _self.push_back(*i_first);
        }
      template<class T>
        bp::object pop_last( std::vector<T> &_self )
          { return pop(_self, _self.size()-1); } 
      template<class T>
        bp::object pop( std::vector<T> &_self, int _i )
        {
          if( _i < 0 ) _i += _self.size();
          if( _i < 0 or _i >= _self.size() )
            python::PyException<error::IndexError>::throw_error("Index out range in pop.");
          bp::object result(_self[_i]);
          _self.erase(_self.begin() + _i);
          return result;
        }
      template<class T>
        void insert( std::vector<T> &_self, int _i, T const &_t )
        {
          if( _i < 0 ) _i += _self.size();
          if( _i < 0 or _i >= _self.size() )
            python::PyException<error::IndexError>::throw_error("Index out range in insert.");
          _self.insert(_self.begin()+_i, _t);
        }
 
      template<class T>
        int index(std::vector<T> &_self, T const &_t )
        {
          typename std::vector<T>::const_iterator i_found = 
             std::find(_self.begin(), _self.end(), _t);
          if(i_found == _self.end()) 
            python::PyException<error::KeyError>::throw_error("Object not in list");
          return i_found - _self.begin();
        }
      template<class T>
        void sort(std::vector<T> &_self)
          { std::sort(_self.begin(), _self.end()); }
      template<class T>
        void reverse(std::vector<T> &_self) 
          { std::reverse(_self.begin(), _self.end()); }
      template<class T>
        int count(std::vector<T> const &_self, T const &_t) 
          { return std::count(_self.begin(), _self.end(), _t); }
      template<class T>
        void remove(std::vector<T> &_self, T const &_t) 
        {
          typename std::vector<T>::iterator i_found = 
             std::find(_self.begin(), _self.end(), _t);
          if(i_found == _self.end()) 
            python::PyException<error::KeyError>::throw_error("Object not in list");
          _self.erase(i_found);
        }
      template<class T>
        bool equal(std::vector<T> const &_a, bp::object _b)
        {
          if(not PyObject_HasAttrString(_b.ptr(), "__iter__"))
            python::PyException<error::AttributeError>::throw_error("Object does not implement __iter__.");
          bp::stl_input_iterator<T> i_first(_b);
          bp::stl_input_iterator<T> const i_end;
          typename std::vector<T>::const_iterator i_in = _a.begin();
          typename std::vector<T>::const_iterator const i_in_end = _a.end();
          for(; i_first != i_end and i_in != i_in_end; ++i_first, ++i_in) 
            if(*i_in != *i_first) return false;
          return i_in == i_in_end and i_first == i_end;
        }
    }

    template<class T>
      bp::class_< std::vector<T> > expose_vector( const std::string &_name,
                                                  const std::string &_docstring )
      {
        return bp::class_< std::vector<T> >( _name.c_str(), _docstring.c_str() )
           .def( "__init__", bp::make_constructor(details::default_constructor<T>) )
           .def( "__init__", bp::make_constructor(details::copy_constructor<T>) )
           .def( "__init__", bp::make_constructor(details::constructor<T>) )
           .def( bp::vector_indexing_suite< std::vector<T>, false >() )
           .def( "__str__", &details::__str__<T> )
           .def( "append", &std::vector<T> :: push_back)
           .def( "extend", &details::extend<T> ) 
           .def( "pop", &details::pop_last<T> ) 
           .def( "pop", &details::pop<T> ) 
           .def( "insert", &details::insert<T> ) 
           .def( "remove", &details::remove<T> ) 
           .def( "index", &details::index<T> ) 
           .def( "reverse", &details::reverse<T> ) 
           .def( "sort", &details::sort<T> ) 
           .def( "count", &details::count<T> ) 
           .def( "__eq__", &details::equal<T> ) 
           .def( "clear", &std::vector<T> :: clear );
      }

  }
} // namespace LaDa
#endif
