//
//  Version: $Id$
//
#include<boost/numeric/ublas/vector_proxy.hpp>
#include<boost/numeric/ublas/matrix_proxy.hpp>
#include<boost/numeric/ublas/operation.hpp>

namespace CE
{
  namespace details 
  {
    template< class T_MATRIX, class T_VECTOR1, class T_VECTOR2, class T_MAPPING > 
      void rank_vector( const T_MATRIX &_mat, const T_VECTOR1 &_in, T_VECTOR2 &_out )
      {
        namespace bblas = boost::numeric::ublas;
        __ASSERT( _in.size() * T_MAPPING::D != _mat.size2(),
                  "Input vector and matrix of incompatible size.\n" )
        __ASSERT( _out.size() != _mat.size1(),
                  "Output vector and matrix of incompatible size.\n" )
        typename T_VECTOR2 :: iterator i_out = _out.begin();
        typename T_VECTOR2 :: iterator i_out_end = _out.end();
        for(size_t i(0); i_out != i_out_end; ++i_out, ++i )
        {
          typedef bblas::matrix_row< T_MATRIX > t_Row;
          t_Row row( _mat(i) );
          *i_out += conf_coef<t_Row, T_VECTOR1, T_MAPPING>( row, _in );
        }
      }
 
    template< class T_VECTOR1, class T_VECTOR2, class T_MAPPING > 
      typename T_VECTOR1::value_type conf_coef( const T_VECTOR1 &_coef,
                                                const T_VECTOR2 &_conf )
      {
        namespace bblas = boost::numeric::ublas;
        __ASSERT( _coef.size() != _conf.size() * T_MAPPING::D,
                  "Coef vector and Conf vector of incompatible size.\n" )
        typename T_VECTOR1 :: value_type result(1);
        typename T_VECTOR2 :: const_iterator i_conf = _conf.begin();
        typename T_VECTOR2 :: const_iterator i_conf_end = _conf.end();
        typename T_VECTOR1 :: const_iterator i_coef = _coef.end();
        for(; i_conf != i_conf_end; ++i_conf, i_coef += T_MAPPING::D )
          T_MAPPING::apply(  *i_conf, i_coef, result );
      }

    template< class T_MATRIX1, class T_MATRIX2, class T_VECTOR, class T_MAPPING > 
      void allconfs_rank_vector( const T_MATRIX1 &_mat, const T_MATRIX2 &_in,
                                 T_VECTOR &_out )
      {
        namespace bblas = boost::numeric::ublas;
        for( size_t i(0); i < _in.size2(); ++i )
        {
          typedef bblas::matrix_column< T_MATRIX2 > t_Column;
          t_Column column( _in( i ) );
          rank_vector<T_MATRIX1, t_Column, T_VECTOR, T_MAPPING>( _mat, column, _out );
        }
      }

    template< class T_MATIN, class T_MATOUT, class T_VEC, class T_MAPPING > 
      void rank_dim_matrix( const T_MATIN &_in, const T_VEC &_vec, T_MATOUT &_out )
      {
        namespace bblas = boost::numeric::ublas;
        __ASSERT( _vec.size() * T_MAPPING::D != _in.size2(),
                  "Incompatible size.\n" )
        __ASSERT( _out.size1() != _vec.size(),
                  "Incompatible size.\n" )
        typename T_VEC :: iterator i_out = _out.begin();
        typename T_VEC :: iterator i_out_end = _out.end();
        for(size_t i(0); i_out != i_out_end; ++i_out, ++i )
        {
          typedef bblas::matrix_row< T_MATIN > t_RowIn;
          typedef bblas::matrix_row< T_MATOUT > t_RowOut;
          t_RowIn rowin( _in(i) );
          t_RowOut rowout( _out(i) );
          conf_coef_vector<t_RowIn, T_VEC,
                           t_RowOut, T_MAPPING>( rowin, _vec, rowout );
        }
      }

    template< class T_VECTOR1, class T_VECTOR2, class T_VECTOR3, class T_MAPPING > 
      typename T_VECTOR1::value_type conf_coef_vector( const T_VECTOR1 &_coef,
                                                       const T_VECTOR2 &_conf,
                                                       const T_VECTOR3 &_out )
      {
        namespace bblas = boost::numeric::ublas;
        __ASSERT( _coef.size() != _conf.size() * T_MAPPING::D,
                  "Coef vector and Conf vector of incompatible size.\n" )
        __ASSERT( _coef.size() != _conf.size(),
                  "Coef vector and Conf vector of incompatible size.\n" )
        typename T_VECTOR3 :: value_type result(0);
        typename T_VECTOR2 :: const_iterator i_conf = _conf.begin();
        typename T_VECTOR2 :: const_iterator i_conf_end = _conf.end();
        typename T_VECTOR1 :: const_iterator i_coef = _coef.end();
        typename T_VECTOR3 :: const_iterator i_out = _out.end();
        for(; i_conf != i_conf_end; ++i_conf, ++i_out, i_coef += T_MAPPING::D )
        {
          typename T_VECTOR3 :: value_type result(1);
          T_MAPPING::apply(  *i_conf, i_coef, result );
          *i_out += result;
        }
      }
  } // end of details namespace.
} // end of CE namespace.
