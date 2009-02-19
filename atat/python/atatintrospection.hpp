//
//  Version: $Id$
//

#ifndef _LADA_PYTHON_ATAT_INTROSPECTION_HPP_
#define _LADA_PYTHON_ATAT_INTROSPECTION_HPP_

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

namespace LaDa
{
  namespace Python
  {
    template< class T_VECTOR > struct vector_introspection;
    template<> struct vector_introspection< atat::rVector3d >
    {
      typedef atat::rVector3d t_Vector;
      typedef types::t_real type;
      const static size_t dim = 3;
      typedef boost::python::init< type, type, type > init;
      struct pickle : boost::python::pickle_suite
      {
        static boost::python::tuple getinitargs( t_Vector const& _w) 
          { return boost::python::make_tuple( _w[0], _w[1], _w[2] ); }
      };

    };
    template<> struct vector_introspection< atat::iVector3d >
    {
      typedef atat::iVector3d t_Vector;
      typedef types::t_int type;
      const static size_t dim = 3;
      typedef boost::python::init< type, type, type > init;
      struct pickle : boost::python::pickle_suite
      {
        static boost::python::tuple getinitargs( t_Vector const& _w) 
          { return boost::python::make_tuple( _w[0], _w[1], _w[2] ); }
      };

    };
    template<> struct vector_introspection< atat::rVector2d >
    {
      typedef atat::rVector2d t_Vector;
      typedef types::t_real type;
      const static size_t dim = 2;
      typedef boost::python::init< type, type > init;
      struct pickle : boost::python::pickle_suite
      {
        static boost::python::tuple getinitargs( t_Vector const& _w) 
          { return boost::python::make_tuple( _w[0], _w[1] ); }
      };
    };
    template<> struct vector_introspection< atat::iVector2d >
    {
      typedef atat::iVector2d t_Vector;
      typedef types::t_int type;
      const static size_t dim = 2;
      typedef boost::python::init< type, type > init;
      struct pickle : boost::python::pickle_suite
      {
        static boost::python::tuple getinitargs( t_Vector const& _w) 
          { return boost::python::make_tuple( _w[0], _w[1] ); }
      };
    };

    template< class T_VECTOR > struct matrix_introspection;
    template<> struct matrix_introspection< atat::rMatrix3d >
    {
      typedef atat::rVector3d t_Vector;
      typedef atat::rMatrix3d t_Matrix;
      typedef types::t_real type;
      const static size_t dim = 3;
      struct pickle : boost::python::pickle_suite
      {
        static boost::python::tuple getinitargs( t_Matrix const& _w)  
        {
          return boost::python::tuple();
        }
        static boost::python::tuple getstate(const t_Matrix& _mat)
        {
          return boost::python::make_tuple( _mat(0,0), _mat(0,1), _mat(0,2),
                                            _mat(1,0), _mat(1,1), _mat(1,2),
                                            _mat(2,0), _mat(2,1), _mat(2,2) );
        }
        static void setstate( t_Matrix& _mat, boost::python::tuple state)
        {
          namespace bp = boost::python;
          if( bp::len( state ) != 9 )
          {
            PyErr_SetObject(PyExc_ValueError,
                            ("expected 9-item tuple in call to __setstate__; got %s"
                             % state).ptr()
                );
            bp::throw_error_already_set();
          }
          size_t u(0);
          for( size_t i(0); i < dim; ++i )
            for( size_t j(0); j < dim; ++j, ++u )
              _mat(i,j) = bp::extract< type >( state[u] );
        }
      };

    };
  }
}

#endif 
