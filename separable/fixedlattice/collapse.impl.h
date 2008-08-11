//
//  Version: $Id$
//
#include<boost/numeric/ublas/vector_proxy.hpp>
#include<boost/numeric/ublas/matrix_proxy.hpp>
#include<boost/numeric/ublas/operation.hpp>

namespace CE
{
# if defined( COLHEAD ) || defined(INCOLLAPSE)
#   error "Macros with same names."
# endif
# define COLHEAD \
    Collapse<T_SEPARABLES, T_MAPPING, T_NORMALIZATION> 
#  define INCOLLAPSE( var ) \
     template< class T_SEPARABLES, class T_MAPPING, class T_NORMALIZATION > \
       var COLHEAD

  INCOLLAPSE(void) :: create_A_n_b( t_Matrix &_A, t_Vector &_b )
  {
    namespace bblas = boost::numeric::ublas;
    __ASSERT( not separables, "Function pointer not set.\n" )
    __ASSERT( dim >= configurations_.size1(), "Inconsistent sizes.\n" )
    // Loop over inequivalent configurations.
    t_Vector X( separables->coefs.size1() * t_Separables::t_Mapping::D );
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

  INCOLLAPSE(void) :: create_X( size_t _i, t_Vector &_out )
  {
    namespace bblas = boost::numeric::ublas;
    __ASSERT(   separables->coefficients.size1()
              * t_Separables::t_Mapping::D == _out.size(),
              "Incompatible sizes.\n" )
    typedef const bblas::matrix_range< t_Matrix > t_Range;
    const bblas::range range( mapping.range( _i ) );
    t_Range confs_range( configurations_, range,
                         bblas::range( 0, configurations_.size1() ) );
    const bblas::matrix_row< t_Range > conf_rows( confs_range, dim );
    typename t_Matrix :: const_iterator i_conf = configurations_.begin();
    typename t_Matrix :: const_iterator i_conf_end = configurations_.end();
    for(size_t c(0); i_conf != i_conf_end; ++i_conf, ++c )
      for(size_t r(0); r < separables->coefficients.size1(); ++r )
      {
        const size_t D( t_Normalization :: D );
        typename t_Vector::value_type scalar( mapping.eweight(_i,c) );
        scalar *= factor( separables->coefficients, range.start() + c, r ); 
        bblas::vector_range< t_Vector > vecr( _out,
                                             bblas::range( r * D, (r+1) * D ) );
        t_Separables::t_Mapping::add_tovec( *i_conf, vecr, scalar );
      }
  }

  INCOLLAPSE( typename COLHEAD :: t_Vector :: value_type )
    :: factor( size_t _kv, size_t _r )
    {
      namespace bblas = boost::numeric::ublas;
      typedef typename t_Separables :: t_Mapping t_SepMap;

      typename t_Vector :: value_type result(1);
      bblas::matrix_row< t_Matrix > row( separables->coefficients, _r );
      t_SepMap :: apply( configurations_( dim, _kv),
                         row.begin() + dim * t_SepMap::D,
                         result );
      if( not Fuzzy::is_zero( result ) ) 
        return scales(_r, _kv) / result * separables->norms[_r];

      // we should never have to be here...
      typedef bblas::matrix_column< t_Matrix > t_Column;
      const t_Column coefs( separables->coefficients, _r );
      const t_Column config( configurations_, _kv );


      result = typename t_Vector::value_type( 1 );
      typename t_Column :: const_iterator i_conf = config.begin();
      typename t_Column :: const_iterator i_conf_end = config.end();
      typename t_Column :: const_iterator i_coef = coefs.end();
      for( size_t d(0); i_conf != i_conf_end; 
           ++d, ++i_conf, i_coef += t_SepMap::D )
        if( d != dim ) t_SepMap::apply( *i_conf, i_coef, result );
      return result * separables->norms[_r];
    }

  INCOLLAPSE( void ) :: update_all()
  {
    namespace bblas = boost::numeric::ublas;
    __ASSERT( scales.size1() == configurations_.size1(),
              "Inconsistent sizes.\n" )
    // First, normalizes coefficients.
    for( size_t r(0); r < separables->coefficients.size1(); ++r )
    {
      typedef bblas::matrix_row< t_Matrix > t_Row;
      t_Row row( separables->coefficients, r );
      typename t_Row :: iterator i_coef = row.begin();
      typename t_Row :: iterator i_coef_end = row.end();
      for(; i_coef != i_coef_end; i_coef += t_Normalization :: D )
        t_Normalization :: apply( i_coef, separables->norms[r] );
    }
    // Then, updates scales. 
    typedef typename t_Matrix :: const_iterator2 :: value_type t_Column;
    typedef typename t_Separables :: t_Mapping t_SepMap;
    typename t_Matrix :: const_iterator2 i_scales = scales.begin2();
    typename t_Matrix :: const_iterator2 i_conf = configurations_.begin2();
    typename t_Matrix :: const_iterator2 i_conf_end = configurations_.end2();
    for(; i_conf != i_conf_end; ++i_conf, ++i_scales )
      details::rank_vector< t_Matrix, t_Column, t_Column, t_SepMap>
                          ( separables->coefficients, *i_conf, *i_scales );
  }

  INCOLLAPSE( template< class T_STRUCTURES> void  )
    :: init( const T_STRUCTURES& _strs, const PosToConfs &_postoconfs )
    {
      namespace bl = boost::lambda;
      __ASSERT( not Crystal::Structure::lattice, 
                "Crystal::Structure::lattice has not been set.\n" )
      // creates configurations.
      typedef std::vector< PosToConfs::t_Configurations > t_Confs; 
      t_Confs confs;
      typename T_STRUCTURES :: const_iterator i_str = _strs.begin();
      typename T_STRUCTURES :: const_iterator i_str_end = _strs.end();
      size_t nbconfs(0);
      for(; i_str != i_str_end; ++i_str )
      {
        PosToConfs::t_Configurations strconf;
        _postoconfs( *i_str, strconf );
        nbconfs += strconf.size();
        confs.push_back( strconf ); 
      }
 
      // translates to matrix.
      configurations_.resize( _postoconfs.positions.size(), nbconfs );
      t_Confs :: const_iterator i_confs = confs.begin();
      t_Confs :: const_iterator i_confs_end = confs.end();
      for(size_t j(0); i_confs != i_confs_end; ++i_confs )
      {
        PosToConfs :: t_Configurations :: const_iterator i_conf = i_confs->begin();
        PosToConfs :: t_Configurations :: const_iterator i_conf_end = i_confs->end();
        for(; i_conf != i_conf_end; ++i_conf, ++j )
        {
          __ASSERT( j != nbconfs, "Inconsistent sizes" );
          for( size_t i(0); i < _postoconfs.positions.size(); ++i )
            configurations_(i,j) = i_conf->first[i] ? 0: 1;
        }
      }
 
      // initializes mapping.
      mapping.init( _strs, confs );
    }

  INCOLLAPSE( opt::ErrorTuple ) :: evaluate()
    {
      opt::ErrorTuple error(0);
      typename t_Matrix :: const_iterator2 i_column = configurations.begin2();
      typename t_Matrix :: const_iterator2 i_column_end = configurations.end2();
      for(size_t n(0); n < mapping.size(); ++n )
      {
        size_t nequiv = mapping.range(n).size();
        types::t_real intermed(0);
        for( size_t j(0); j < nequiv; ++j, ++i_column )
          intermed += (*separables)( *i_column ) * mapping.eweights(n,j);
        error += ErrorTuple( mapping.targets - intermed, mapping.weights(n) );
      }
      return error;
    }

# undef COLHEAD
# undef INCOLLAPSE
} // end of CE namespace.
