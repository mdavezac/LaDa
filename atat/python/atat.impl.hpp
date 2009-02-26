//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <iostream>
#include <algorithm>
#include <stdexcept>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/python.hpp>
#include <boost/python/object.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/list.hpp>
#include <boost/lexical_cast.hpp>

#include <opt/types.h>
#include <opt/debug.h>

#define STREAM_VECTOR
#include "../fxvector.h"
#include "../vectmac.h"


#include "atatintrospection.hpp"
#include "atatoperators.hpp"

namespace LaDa
{
  namespace Python
  {
    namespace details 
    {
      template< class T_VECTOR > 
        typename vector_introspection< T_VECTOR > :: type 
          norm2( const T_VECTOR &_vector ) { return LaDa::atat::norm2( _vector ); }
 
      template< class T_VECTOR >
        typename vector_introspection<T_VECTOR>::type
          getvecitem( const T_VECTOR& _vec, types::t_int _i )
          {
            typedef typename vector_introspection< T_VECTOR > :: type type;
            const types::t_int dim(  vector_introspection< T_VECTOR > :: dim );
            if( _i > dim or _i < -dim )
            {
              std::ostringstream sstr;
              throw std::out_of_range( "atat vector" );
            }
            return _vec.x[ size_t( _i < 0 ? dim + _i: _i ) ];
          }
      template< class T_VECTOR >
        void setvecitem( T_VECTOR& _vec, types::t_int _i,
                         typename vector_introspection< T_VECTOR > :: type _a )
        {
          typedef typename vector_introspection< T_VECTOR > :: type type;
          const types::t_int dim(  vector_introspection< T_VECTOR > :: dim );
          if( _i > dim or _i < -dim )
          {
            std::ostringstream sstr;
            throw std::out_of_range( "atat vector" );
          }
          _vec.x[ size_t( _i < 0 ? types::t_int(dim) + _i: _i ) ] = _a;
        }
      template< class T_VECTOR >
        size_t getlength( const T_VECTOR& _vec )
          { return vector_introspection<T_VECTOR> :: dim; }
      template< class T_VECTOR >
        std::string print( const T_VECTOR& _vec )
        { 
          std::ostringstream sstr;
          sstr << _vec;
          return sstr.str();
        }
 
      template< class T_MATRIX >
        typename matrix_introspection<T_MATRIX>::type
          getmatitem( const T_MATRIX& _m, const boost::python::tuple &_t )
          {
            namespace bp = boost::python;
            const types::t_int _i = bp::extract< types::t_int >( _t[0] );
            const types::t_int _j = bp::extract< types::t_int >( _t[1] );
            const types::t_int dim(  matrix_introspection< T_MATRIX > :: dim );
            if( _i > dim or _i < -dim or _j > dim or _j < -dim )
            {
              std::ostringstream sstr;
              throw std::out_of_range( "atat matrix" );
            }
            return _m.x[ size_t(_i < 0 ? dim + _i: _i) ][ size_t(_j < 0 ? dim + _j: _j) ];
          }
      template< class T_MATRIX >
        void setmatitem( T_MATRIX& _m, const boost::python::tuple &_t,
                         typename matrix_introspection< T_MATRIX > :: type _a )
        {
          namespace bp = boost::python;
          const types::t_int _i = bp::extract< types::t_int >( _t[0] );
          const types::t_int _j = bp::extract< types::t_int >( _t[1] );
          const types::t_int dim(  matrix_introspection< T_MATRIX > :: dim );
          if( _i > dim or _i < -dim or _j > dim or _j < -dim )
          {
            std::ostringstream sstr;
            throw std::out_of_range( "atat matrix" );
          }
          _m.x[ size_t(_i < 0 ? dim + _i: _i) ][ size_t(_j < 0 ? dim + _j: _j) ] = _a;
        }
      template< class T_MATRIX >
        size_t getmlength( const T_MATRIX& )
          { return matrix_introspection<T_MATRIX> :: dim; }
    } // namespace details

    template< class T_VECTOR >
      T_VECTOR* empty_vector()
      {
        namespace bp = boost::python;
        typedef typename vector_introspection< T_VECTOR > :: type type;
        const size_t dim(  vector_introspection< T_VECTOR > :: dim );
        return new T_VECTOR( 0,0,0 );
      }
    
    template< class T_VECTOR >
      T_VECTOR* copy_vector( const T_VECTOR &_o )
      {
        namespace bp = boost::python;
        typedef typename vector_introspection< T_VECTOR > :: type type;
        const size_t dim(  vector_introspection< T_VECTOR > :: dim );
        return new T_VECTOR( _o );
      }

    template< class T_VECTOR >
      T_VECTOR* make_vector( boost::python::object &_o )
      {
        namespace bp = boost::python;
        typedef typename vector_introspection< T_VECTOR > :: type type;
        const size_t dim(  vector_introspection< T_VECTOR > :: dim );
        T_VECTOR *result = new T_VECTOR;
        try
        {
          const size_t size( bp::len( _o ) );
          __DOASSERT( size != dim and size != 1, "Incorrect size.\n" )
          for( size_t i=0; i < dim; ++i )
            (*result)(i) = bp::extract< type >( _o[i] );
          return result;
        }
        catch( ... )
        {
          delete result;
          __DOASSERT( true, "" ) 
        }
      }

    template< class T_MATRIX >
      T_MATRIX* empty_matrix()
      {
        namespace bp = boost::python;
        typedef typename matrix_introspection< T_MATRIX > :: type type;
        const size_t dim(  matrix_introspection< T_MATRIX > :: dim );
        T_MATRIX *result = new T_MATRIX; result->zero();
        return result; 
      }
    template< class T_MATRIX >
      T_MATRIX* copy_matrix( const T_MATRIX &_o )
      {
        namespace bp = boost::python;
        typedef typename matrix_introspection< T_MATRIX > :: type type;
        const size_t dim(  matrix_introspection< T_MATRIX > :: dim );
        return new T_MATRIX(_o);
      }
    template< class T_MATRIX >
      T_MATRIX* make_matrix( boost::python::object &_o )
      {
        namespace bp = boost::python;
        typedef typename matrix_introspection< T_MATRIX > :: type type;
        const size_t dim(  matrix_introspection< T_MATRIX > :: dim );
        T_MATRIX *result = new T_MATRIX; result->zero();
        try
        {
          const size_t size( bp::len( _o ) );
          __DOASSERT( size != dim and size != 1, "Incorrect size.\n" )
          for( size_t i=0; i < dim; ++i )
          {
            const boost::python::object object = bp::extract<bp::object>( _o[i] );
            const size_t inner = bp::len(_o[i]); 
            __DOASSERT( inner != dim, "Incorrect size.\n" )
            for( size_t j=0; j < dim; ++j )
              (*result)(i,j) = bp::extract< type >( object[j] );
          }
          return result;
        }
        catch( ... )
        {
          delete result;
          __DOASSERT( true, "" ) 
        }
      }

    template< class T_VECTOR >
      void expose_atatvector( const std::string &_name, const std::string &_docstring )
      {
        namespace bp = boost::python;
        typedef typename vector_introspection< T_VECTOR > :: type type;
        typedef typename vector_introspection< T_VECTOR > :: pickle pickle;
        typedef typename vector_introspection< T_VECTOR > :: init init;
        const size_t dim(  vector_introspection< T_VECTOR > :: dim );
        // expose atat vector.
        bp::class_< T_VECTOR>( _name.c_str(), _docstring.c_str() )
            .def( "__init__", bp::make_constructor( make_vector<T_VECTOR> ) )
            .def( "__init__", bp::make_constructor( copy_vector<T_VECTOR> ) )
            .def( "__init__", bp::make_constructor( empty_vector<T_VECTOR> ) )
            .def( bp::self + bp::other< T_VECTOR >() ) 
            .def( bp::other< T_VECTOR >() + bp::self ) 
            .def( bp::self - bp::other< T_VECTOR >() ) 
            .def( bp::other< T_VECTOR >() - bp::self ) 
            .def( bp::other< T_VECTOR >() * bp::self ) 
            .def( bp::self * bp::other< T_VECTOR >() )
            .def( type() * bp::self )
            .def( bp::self * type() )
            .def( bp::self / type() )
            .def( bp::self == bp::other<T_VECTOR>() )
            .def( bp::self != bp::other<T_VECTOR>() )
            .def( "__getitem__", &details::getvecitem<T_VECTOR> )
            .def( "__setitem__", &details::setvecitem<T_VECTOR> ) 
            .def( "__len__", &details::getlength<T_VECTOR> ) 
            .def( "__str__", &details::print<T_VECTOR> )
            .def_pickle( pickle() );

        bp::def( _name.c_str(), &make_vector< T_VECTOR >, 
                 bp::return_value_policy<bp::manage_new_object>(),
                 _docstring.c_str() );
        bp::def( "norm2", &details::norm2<T_VECTOR>,
                 bp::arg("vec"),
                 ("Returns squared euclidian-norm of an " + _name + " object.").c_str() );
      }

    template< class T_MATRIX >
      void expose_atatmatrix( const std::string &_name, const std::string &_docstring )
      {
        namespace bp = boost::python;
        typedef typename matrix_introspection< T_MATRIX > :: type type;
        typedef typename matrix_introspection< T_MATRIX > :: t_Vector t_Vector;
        typedef typename matrix_introspection< T_MATRIX > :: pickle pickle;
        const size_t dim(  matrix_introspection< T_MATRIX > :: dim );
        // expose atat vector.
        bp::class_< T_MATRIX>( _name.c_str(), _docstring.c_str() )
            .def( "__init__", bp::make_constructor( make_matrix<T_MATRIX> ) )
            .def( "__init__", bp::make_constructor( copy_matrix<T_MATRIX> ) )
            .def( "__init__", bp::make_constructor( empty_matrix<T_MATRIX> ) )
            .def( bp::self + bp::other< T_MATRIX >() ) 
            .def( bp::other< T_MATRIX >() + bp::self ) 
            .def( bp::self - bp::other< T_MATRIX >() ) 
            .def( bp::other< T_MATRIX >() - bp::self ) 
            .def( bp::other< T_MATRIX >() * bp::self ) 
            .def( bp::self * bp::other< T_MATRIX >() )
            .def( bp::self * bp::other< t_Vector >() )
            .def( bp::other< t_Vector >() * bp::self )
            .def( type() * bp::self )
            .def( bp::self * type() )
            .def( bp::self / type() )
            .def( bp::self != bp::other<T_MATRIX>() )
            .def( bp::self == bp::other<T_MATRIX>() )
            .def( "__getitem__", &details::getmatitem<T_MATRIX> )
            .def( "__setitem__", &details::setmatitem<T_MATRIX> ) 
            .def( "__len__", &details::getmlength<T_MATRIX> ) 
            .def( "__str__", &details::print<T_MATRIX> )
            .def_pickle( pickle() );
      }

  }
}

