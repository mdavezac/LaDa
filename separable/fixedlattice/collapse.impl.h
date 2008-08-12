//
//  Version: $Id$
//
#include<boost/numeric/ublas/vector_proxy.hpp>
#include<boost/numeric/ublas/matrix_proxy.hpp>
#include<boost/numeric/ublas/operation.hpp>
//#include<boost/lambda/if.hpp>

namespace CE
{
# if defined( COLHEAD ) || defined(INCOLLAPSE)
#   error "Macros with same names."
# endif
# define COLHEAD \
    Collapse<T_SEPARABLES, T_MAPPING> 
#  define INCOLLAPSE( var ) \
     template< class T_SEPARABLES, class T_MAPPING >  var COLHEAD

  INCOLLAPSE(void) :: create_A_n_b( t_Matrix &_A, t_Vector &_b )
  {
    namespace bblas = boost::numeric::ublas;
    __ASSERT( not separables, "Function pointer not set.\n" )
    __ASSERT( dim >= configurations_.size1(), "Inconsistent sizes.\n" )
    // Loop over inequivalent configurations.
    t_Vector X( separables->coefficients.size1() );
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
    __ASSERT(   separables->coefficients.size1() != _out.size(),
              "Incompatible sizes.\n" )
    // Create matrix range including only equivalent configurations.
    typedef const bblas::matrix_range< t_Matrix > t_Range;
    const bblas::range equivrange( mapping.range( _i ) );
    const bblas::range dimrange( dim, dim+1 );
    t_Range equivconfs( configurations_, dimrange, equivrange );

    typename t_Range :: const_iterator2 i_conf = equivconfs.begin2();
    typename t_Range :: const_iterator2 i_conf_end = equivconfs.end2();
    for(size_t c(0); i_conf != i_conf_end; ++i_conf, ++c )
      for(size_t r(0); r < separables->ranks(); ++r )
      {
        const size_t D( t_Separables::t_Mapping :: D );
        typename t_Vector::value_type scalar( mapping.eweight(_i,c) );
        scalar *= factor( equivrange.start() + c, r ); 
        const bblas::range range( r * D, (r+1) * D );
        bblas::vector_range< t_Vector > vecr( _out, range );
        t_Separables::t_Mapping::add_tovec( *i_conf, vecr, scalar );
      }
  }

  INCOLLAPSE( typename COLHEAD :: t_Vector :: value_type )
    :: factor( size_t _kv, size_t _r )
    {
      namespace bl = boost::lambda;
      namespace bblas = boost::numeric::ublas;

      typename t_Vector :: value_type result;
      bblas::matrix_column< t_Matrix > config( configurations_, _kv );
      t_Separables :: t_Policy :: apply_to_dim_n_rank( separables->coefficients, config, 
                                                       result, dim, _r, bl::_1 = bl::_2 );
      if( not Fuzzy::is_zero( result ) ) 
        return scales(_r, _kv) / result * separables->norms[_r];

      // we should never have to be here...
      result = typename t_Vector :: value_type(1);
      size_t d(0);
      t_Separables :: t_Policy :: apply_to_rank
      ( 
        separables->coefficients, config, result, _r,
        (
          bl::if_then( bl::var(d) != bl::constant(dim), bl::_1 *= bl::_2 ),
          ++bl::var(d)
        )
      );
      return result * separables->norms[_r];
    }

  INCOLLAPSE( void ) :: update_all()
  {
    namespace bl = boost::lambda;
    namespace bblas = boost::numeric::ublas;
    __ASSERT( scales.size1() == configurations_.size1(),
              "Inconsistent sizes.\n" )
    // First, normalizes coefficients.
    separables->normalize();

    // Then, updates scales. 
    for( size_t i(0); i < configurations_.size2(); ++i )
    {
      bblas::matrix_column<t_Matrix> conf( configurations_, i );
      bblas::matrix_column<t_Matrix> scaling( scales, i );
      std::fill( scaling.begin(), scaling.end(), typename t_Matrix::value_type(1) );
      t_Separables::t_Policy::rank_vector
      ( 
        separables->coefficients, conf, scaling,
        bl::_1 *= bl::_2
      );
    }
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
          __ASSERT( j == nbconfs, "Inconsistent sizes" );
          for( size_t i(0); i < _postoconfs.positions.size(); ++i )
            configurations_(i,j) = i_conf->first[i] ? 0: 1;
        }
      }
 
      // initializes mapping.
      mapping.init( _strs, confs );
    }

  INCOLLAPSE( opt::ErrorTuple ) :: evaluate()
  {
    namespace bblas = boost::numeric::ublas;
    opt::ErrorTuple error(0);
    for(size_t n(0); n < mapping.size(); ++n )
    {
      bblas::range range( mapping.range(n) );
      types::t_real intermed(0);
      for( bblas::range::const_iterator j( range.begin() ); j != range.end(); ++j )
      {
        const bblas::matrix_column<t_Matrix> config( configurations_, *j );
        intermed += (*separables)( config ) * mapping.eweight(n,*j);
      }
      error += opt::ErrorTuple( mapping.target(n) - intermed, mapping.weight(n) );
    }
    return error;
  }

  INCOLLAPSE( void ) :: init( t_Separables& _sep )
  {
    separables = &_sep; 
    scales.resize( separables->coefficients.size1() / t_Separables::t_Mapping::D,
                   configurations_.size2() );
  }
# undef COLHEAD
# undef INCOLLAPSE
} // end of CE namespace.
