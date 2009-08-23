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
#include <boost/python/class.hpp>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <opt/debug.h>

namespace LaDa
{
  namespace Python
  {
    namespace details
    {
      template< class T_TYPE > std::vector<T_TYPE>* default_constructor();
      template< class T_TYPE > 
        std::vector<T_TYPE>* copy_constructor( const std::vector<T_TYPE>& _ob );
      template< class T_TYPE >
        std::vector<T_TYPE>* object_constructor( const boost::python::list& _ob );
      template< class T_TYPE >
        std::string print( const std::vector<T_TYPE>& _ob );
    }

    template< class T_TYPE >
      void expose_vector( const std::string &_name, const std::string &_docstring )
      {
        namespace bp = boost::python;
        bp::class_< std::vector< T_TYPE > >( _name.c_str(), _docstring.c_str() )
          .def( "__init__", bp::make_constructor( details::default_constructor< T_TYPE > ) )
          .def( "__init__", bp::make_constructor( details::copy_constructor< T_TYPE > ) )
          .def( "__init__", bp::make_constructor( details::object_constructor< T_TYPE > ) )
          .def( bp::vector_indexing_suite< std::vector<T_TYPE> >() )
#         ifndef LADA_PYTHON_STD_VECTOR_NOPRINT
            .def( "__str__", &details::print<T_TYPE> )
#         else
#           undef rADA_PYTHON_STD_VECTOR_NOPRINT
#         endif
          .def( "clear", &std::vector<T_TYPE> :: clear );
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
            sstr << *i_var << " ";
          return sstr.str();
        }
    }
  }
} // namespace LaDa
#endif
