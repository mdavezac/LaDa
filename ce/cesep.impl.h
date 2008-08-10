//
//  Version: $Id$
//
#include<boost/blas/numeric/vector_proxies.hpp>
#include<boost/blas/numeric/matrix_proxies.hpp>
#include<boost/blas/numeric/operations.hpp>

namespace CE
{
  template< class T_MATRIX, class T_SEPS >
  types::t_real sum_over_confs( const T_SEPS &_seps, const T_MATRIX &_confs )
  {
    t_Vector intermed( coefficients.size1() );
    std::fill( intermed.begin(), intermed.end(), 0e0 );
    details::allconfs_rank_vector< typename T_SEPS :: t_Matrix,
                                   T_MATRIX, 
                                   typename T_SEPS :: t_Vector,
                                   typename T_SEPS :: t_Mapping >
                                 ( _seps.coefficients, _confs, intermed );
    return bblas::inner_prod( intermed, _seps.norms );
  }

  template< class T_SEPARABLES, T_MAPPING >
    void Collapse<T_SEPARABLES, T_MAPPING>
      :: create_A_n_b( t_Matrix &_A, t_Vector &_b, const t_Matrix &_coef )
      {
        __ASSERT( dim >= configurations.size1(),
                  "Inconsistent sizes.\n" )
        // Loop over inequivalent configurations.
        t_Vector X( nb_ranks * t_Separable::t_Mapping::D );
        for( size_t i(0); i < mapping.size(); ++i )
        {
          // allows leave-one-out, or leave-many-out.
          if( mapping.do_skip(i) ) continue;
     
          // create the X vector.
          std::fill( X.begin(), X.end(), typename t_Vector::value_type(0) );
          create_X( i, X );
     
          _A += mapping.weight(i) * bblas::outer_prod( X, X ); 
          _b += mapping.weight(i) * mapping.target(i) * X;
        }
      }

  template< class T_SEPARABLES, T_MAPPING >
    void Collapse<T_SEPARABLES, T_MAPPING> :: create_X( size_t _i,
                                                        t_Vector &_out,
                                                        const t_Matrix &_coefs  )
    {
      namespace bblas = boost::numeric::ublas;
      __ASSERT( _nbranks * t_Separables::t_Mapping::D == _out.size(),
                "Incompatible sizes.\n" )
      typedef const bblas::matrix_range< t_Matrix > t_Range;
      const bblas::range range( mapping.range( _i ) );
      t_Range confs_range( configurations, range,
                           bblas::range( 0, configurations.size1() ) );
      const matrix_row< t_Range > conf_rows( confs_range, dim );
      typename T_CONFS :: const_iterator i_conf = confs.begin();
      typename T_CONFS :: const_iterator i_conf_end = confs.end();
      for(size_t c(0); i_conf != i_conf_end; ++i_conf, ++c )
        for(size_t r(0); r < _coefs.size1(); ++r )
        {
          typename t_Vector::value_type scalar( mapping.eweight(_i,c) );
          scalar *= factor( _coef, range.start() + c, r ); 
          bblas::vector_range< t_Vector > vecr( _out, range( r * D, (r+1) * D ) );
          t_Separables::t_Mapping::add_tovec( *i_conf, vecr, scalar );
        }
    }

  template< class T_SEPARABLES, T_MAPPING >
    typename T_SEPARABLES::t_Vector::value_type
      Collapse<T_SEPARABLES, T_MAPPING>  
        :: factor( t_Matrix &_coef, size_t _kv, size_t _r )
        {
          namespace bblas = boost::numeric::ublas;
          typedef typename t_Separables :: t_Mapping t_SepMap;

          typename t_Vector :: value_type result(1);
          bblas::matrix_row< t_Matrix > row( _coef, _r );
          t_SepMap :: apply( configurations( dim, _kv),
                             row.begin() + dim * t_SepMap::D,
                             result );
          if( not Fuzzy::is_zero( result ) ) 
            return scaling(_r, _kv) / result;

          // we should never have to be here...
          typedef bblas::matrix_column< t_Matrix > t_Column;
          const t_Column coefs( _coefs, _r );
          const t_Column config( configurations, _kv );


          result = typename t_Vector::value_type( 1 );
          t_Column :: const_iterator i_conf = config.begin();
          t_Column :: const_iterator i_conf_end = config.end();
          t_Column :: const_iterator i_coef = coefs.end();
          for( size_t d(0); i_conf != i_conf_end; 
               ++d, ++i_conf, i_coef += t_SepMap::D )
            if( d != dim ) t_SepMap::apply( *i_conf, i_coef, result );
          return result;
        }

  template< class T_SEPARABLES, T_MAPPING >
    void Collapse<T_SEPARABLES, T_MAPPING> :: update_all( const t_Matrix &_coefs )
      {
        __ASSERT( scales.size1() == configurations.size1(),
                  "Inconsistent sizes.\n" )
        typedef typename t_Matrix :: const_iterator2 :: value_type t_Column;
        typedef typename t_Separables :: t_Mapping t_SepMap;
        t_Matrix :: const_iterator2 i_scales = scales.begin2();
        t_Matrix :: const_iterator2 i_conf = configurations.begin2();
        t_Matrix :: const_iterator2 i_conf_end = configurations.end2();
        for(; i_conf != i_conf_end; ++i_conf, ++i_scales )
          details::rank_vector< t_Matrix, t_Column, t_Column, t_SepMap>
                              ( _coefs, *i_conf, *i_scales );
      }

  template< class T_SEPARABLES, T_MAPPING > template< class T_STRUCTURES >
    void Collapse<T_SEPARABLES, T_MAPPING> :: init( const T_STRUCTURES& _strs, 
                                                    const std::string& _bdesc )
    {
      namespace bl = boost::lambda;
      __ASSERT( not Crystal::Structure::lattice, 
                "Crystal::Structure::lattice has not been set.\n" )
      // creates configurations.
      PosToConfs postoconfs( *Crystal::Structure::lattice );
      typedef std::vector< PosToConfs::t_Configuations > t_Confs; 
      t_Confs confs;
      postoconf.create_positions( _bdesc );
      typename T_STRUCTURES :: const_iterator i_str = structures.begin();
      typename T_STRUCTURES :: const_iterator i_str_end = structures.end();
      size_t nbconfs(0);
      for(; i_str != i_str_end; ++i_str )
      {
        PosToConfs strconf;
        nbconfs += strconf.size();
        postoconfs( *i_str, strconf );
        confs.push_back( strconf ); 
      }

      // translates to matrix.
      configurations.resize( postoconfs.basis.size(), nbconfs );
      t_Configurations :: const_iterator i_confs = confs.begin();
      t_Configurations :: const_iterator i_confs_end = confs.end();
      for(size_t j(0); i_confs != i_confs_end; ++i_confs )
      {
        PosToConfs :: t_Configurations :: const_iterator i_conf = i_confs->begin();
        PosToConfs :: t_Configurations :: const_iterator i_conf_end = i_confs->end();
        for(; i_conf != i_conf_end; ++i_conf, ++j )
        {
          __ASSERT( j != nbconfs, "Inconsistent sizes" );
          for( size_t i(0); i < postoconfs.basis.size(); ++i )
            configurations(i,j) = i_conf->first[i] ? 0: 1;
        }
      }

      // initializes mapping.
      mapping.init( _strs, confs );
    }

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
        __ASSERT( _coef.size() != conf.size() * T_MAPPING::D,
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
        for( size_t i(0); i < _in.size2(); ++i )
        {
          typedef bblas::matrix_column< T_MATRIX > t_Column;
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
        typename T_VECTOR2 :: iterator i_out = _out.begin();
        typename T_VECTOR2 :: iterator i_out_end = _out.end();
        for(size_t i(0); i_out != i_out_end; ++i_out, ++i )
        {
          typedef bblas::matrix_row< T_MATIN > t_RowIn;
          typedef bblas::matrix_row< T_MATOUT > t_RowOut;
          t_RowIn rowin( _in(i) );
          t_RowOut rowout( _out(i) );
          conf_coef_vector<t_RowIn, T_VEC, t_RowOut, T_MAPPING>( rowin, _vec, rowout )
        }
      }

    template< class T_VECTOR1, class T_VECTOR2, class T_VECTOR3, class T_MAPPING > 
      typename T_VECTOR1::value_type conf_coef_vector( const T_VECTOR1 &_coef,
                                                       const T_VECTOR2 &_conf,
                                                       const T_VECTOR3 &_out )
      {
        namespace bblas = boost::numeric::ublas;
        __ASSERT( _coef.size() != conf.size() * T_MAPPING::D,
                  "Coef vector and Conf vector of incompatible size.\n" )
        __ASSERT( _coef.size() != conf.size(),
                  "Coef vector and Conf vector of incompatible size.\n" )
        typename T_VECTOR :: value_type result(0);
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
