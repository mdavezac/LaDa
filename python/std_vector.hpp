//
//  Version: $Id$
//

#ifndef _LADA_OPT_PYTHON_STD_VECTOR_HPP_
#define _LADA_OPT_PYTHON_STD_VECTOR_HPP_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
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

#include <opt/debug.h>
#include <python/iterator.hpp>

namespace LaDa
{
  namespace Python
  {
    namespace bp = boost::python;
    namespace details
    {
      template< class T_TYPE > std::vector<T_TYPE>* default_constructor();
      template< class T_TYPE > 
        std::vector<T_TYPE>* copy_constructor( const std::vector<T_TYPE>& _ob );
      template< class T_TYPE >
        std::vector<T_TYPE>* object_constructor( const boost::python::list& _ob );
      template< class T_TYPE >
        std::string print( const std::vector<T_TYPE>& _ob );
      template<class T>
        void extend( std::vector<T> & _self, bp::object const &_object )
        {
          size_t const N = bp::len(_object);
          for(size_t i(0); i < N; ++i)
            try { _self.push_back(  bp::extract<T>(_object[i]) ); }
            catch(std::exception &_e)
            {
              PyErr_SetString
              ( 
                PyExc_RuntimeError,
                ("Could not extract object.\n" + std::string(_e.what())).c_str()
              );
              bp::throw_error_already_set();
            }
        }
      template<class T>
        bp::object pop_last( std::vector<T> &_self )
          { return pop(_self, -1); } 
      template<class T>
        bp::object pop( std::vector<T> &_self, int _i )
        {
          if( _i < 0 ) _i += _self.size();
          if( _i < 0 or _i >= _self.size() )
          {
            PyErr_SetString(PyExc_IndexError, "Index out-of-range.\n");
            bp::throw_error_already_set();
            return bp::object();
          }
          typename std::vector<T> :: iterator i_var = _self.begin();
          for(size_t i(0); i < _i; ++i) ++i_var; 
          bp::object result(*i_var);
          _self.erase(i_var);
          return result;
        }
      template<class T>
        void insert( std::vector<T> &_self, T const &_t, size_t _i )
        {
          if( _i < 0 ) _i += _self.size();
          if( _i < 0 or _i >= _self.size() )
          {
            PyErr_SetString(PyExc_IndexError, "Index out-of-range.\n");
            bp::throw_error_already_set();
            return;
          }
          typename std::vector<T> :: iterator i_var = _self.begin();
          for(size_t i(0); i < _i; ++i) ++i_var; 
          _self.insert(i_var, _t);
        }

    }

    template< class T_TYPE >
      bp::class_< std::vector<T_TYPE> > expose_vector( const std::string &_name,
                                                       const std::string &_docstring )
      {
        typedef typename boost::mpl::if_
             < 
               boost::is_class<T_TYPE>,
               bp::return_internal_reference<>,
               bp::return_value_policy<bp::return_by_value>
             > :: type t_policies;
        bp::class_< std::vector<T_TYPE> > result( _name.c_str(), _docstring.c_str() );
        result
          .def( "__init__", bp::make_constructor( details::default_constructor< T_TYPE > ) )
          .def( "__init__", bp::make_constructor( details::copy_constructor< T_TYPE > ) )
          .def( "__init__", bp::make_constructor( details::object_constructor< T_TYPE > ) )
          .def( bp::vector_indexing_suite< std::vector<T_TYPE>, false >() )
#         ifndef LADA_PYTHON_STD_VECTOR_NOPRINT
            .def( "__str__", &details::print<T_TYPE> )
#         else
#           undef rADA_PYTHON_STD_VECTOR_NOPRINT
#         endif
          .def( "append", &std::vector<T_TYPE> :: push_back, bp::with_custodian_and_ward<1,2>() ) 
          .def( "extend", &details::extend<T_TYPE> ) 
          .def( "pop", &details::pop_last<T_TYPE> ) 
          .def( "pop", &details::pop<T_TYPE> ) 
          .def( "insert", &details::insert<T_TYPE>, bp::with_custodian_and_ward<1,2>() ) 
          .def( "clear", &std::vector<T_TYPE> :: clear );
        return result;
      }

    namespace details
    {
      template< class T_TYPE > std::vector<T_TYPE>* default_constructor()
        { return new std::vector<T_TYPE>; }
      template< class T_TYPE > 
        std::vector<T_TYPE>* copy_constructor( const std::vector<T_TYPE>& _ob )
        { return new std::vector<T_TYPE>(_ob); }
      template< class T_TYPE >
        std::vector<T_TYPE>* object_constructor( const boost::python::list& _ob )
        {
          namespace bp = boost::python;
          const size_t n( bp::len( _ob ) );
          std::vector<T_TYPE>* result = new std::vector<T_TYPE>;
          result->reserve( n );
          for( size_t i(0); i < n; ++i )
            result->push_back( bp::extract<T_TYPE>( _ob[i] ) );
          return result;
        }

      template< class T_TYPE >
        std::string print( const std::vector<T_TYPE>& _ob )
        {
          std::ostringstream sstr;
          typename std::vector<T_TYPE>::const_iterator i_var = _ob.begin();
          typename std::vector<T_TYPE>::const_iterator i_var_end = _ob.end();
          for(; i_var != i_var_end; ++i_var )
            sstr << (*i_var) << " ";
          return sstr.str();
        }
    }
  }
} // namespace LaDa
#endif
