//
//  Version: $Id$
//
#include<boost/numeric/ublas/vector_proxy.hpp>
#include<boost/numeric/ublas/matrix_proxy.hpp>
#include<boost/numeric/ublas/operation.hpp>
#include<boost/lambda/bind.hpp>

namespace CE
{
# if defined( COLHEAD ) || defined(INCOLLAPSE) || defined(INCOLLAPSE2)
#   error "Macros with same names."
# endif
# define COLHEAD \
    Collapse<T_TRAITS> 
#  define INCOLLAPSE( var ) \
     template< class T_TRAITS >  var COLHEAD
#  define INCOLLAPSE2( code1, code2 ) \
     template< class T_TRAITS >  code1, code2 COLHEAD

  INCOLLAPSE2(template< class T_MATRIX, class T_VECTOR > void) 
    :: create_A_n_b( T_MATRIX &_A, T_VECTOR &_b )
    {
      namespace bblas = boost::numeric::ublas;
      __ASSERT( not separables_, "Function pointer not set.\n" )
      __ASSERT( dim >= configurations_->size1(), "Inconsistent sizes.\n" )
      // Loop over inequivalent configurations.
      t_Vector X( dof() );
      std::fill( _A.data().begin(), _A.data().end(), 
                 typename t_Matrix::value_type(0) );
      std::fill( _b.data().begin(), _b.data().end(), 
                 typename t_Vector::value_type(0) );
      for( size_t i(0); i < mapping().size(); ++i )
      {
        // allows leave-one-out, or leave-many-out.
        if( mapping().do_skip(i) ) continue;
    
        // create the X vector.
        std::fill( X.begin(), X.end(), typename t_Vector::value_type(0) );
        create_X( i, X );
        std::cout << X << "\n";
        std::cout << coefficients() << "\n";
    
        _A += mapping().weight(i) * bblas::outer_prod( X, X ); 
        _b += mapping().weight(i) * mapping().target(i) * X;
      }
    }

  INCOLLAPSE(template< class T_VECTOR > void)
    :: create_X( size_t _i, T_VECTOR &_out )
    {
      namespace bblas = boost::numeric::ublas;
      __ASSERT( dof() != _out.size(), "Incompatible sizes.\n" )
      // Create matrix range including only equivalent configurations.
      typedef const bblas::matrix_range< t_Configurations > t_Range;
      const bblas::range equivrange( mapping().range( _i ) );
      const bblas::range dimrange( dim, dim+1 );
      t_Range equivconfs( *configurations_, dimrange, equivrange );
 
      typename t_Range :: const_iterator2 i_conf = equivconfs.begin2();
      typename t_Range :: const_iterator2 i_conf_end = equivconfs.end2();
      for(size_t c(0); i_conf != i_conf_end; ++i_conf, ++c )
        for(size_t r(0); r < separables().ranks(); ++r )
        {
          const size_t D( t_Separables::t_Mapping :: D );
          typename t_Vector::value_type scalar( mapping().eweight(_i,c) );
          scalar *=   update_.factor( equivrange.start() + c, r, dim )
                    * separables().norms[r]; 
          const bblas::range range( r * D, (r+1) * D );
          bblas::vector_range< T_VECTOR > vecr( _out, range );
          t_Separables::t_Mapping::add_tovec( *i_conf, vecr, scalar );
        }
    }

  INCOLLAPSE(void) :: update_all()
  {
    update_.init( configurations_ );
    // First, normalizes coefficients.
    separables().normalize();
    // then calls policy.
    update_();
  }
  INCOLLAPSE(void) :: update( types::t_unsigned _d )
  {
    // First, normalizes coefficients.
    separables().normalize();
    // then calls policy.
    update_( _d );
  }
           


  INCOLLAPSE( template< class T_STRUCTURES> void  )
    :: init( const T_STRUCTURES& _strs, const PosToConfs &_postoconfs )
    {
      namespace bblas = boost::numeric::ublas;
      __ASSERT( not Crystal::Structure::lattice, 
                "Crystal::Structure::lattice has not been set.\n" )
      // creates configurations.
      typedef std::vector< PosToConfs::t_Configurations > t_Confs; 
      t_Confs confs;
      size_t nbconfs(0);
      foreach( const typename T_STRUCTURES::value_type& _structure, _strs )
      {
        PosToConfs::t_Configurations strconf;
        _postoconfs( _structure, strconf );
        nbconfs += strconf.size();
        confs.push_back( strconf ); 
      }
 
      // translates to matrix.
      configurations_->resize( _postoconfs.positions.size(), nbconfs );
      size_t i(0);
      foreach( t_Confs :: value_type& confs_per_str, confs )
      {
        foreach( t_Confs :: value_type :: value_type& config, confs_per_str )
        {
          bblas::matrix_column< t_Configurations > column( configurations(), i );
          std::copy( config.first.begin(), config.first.end(), column.begin() );
          ++i;
        }
      }
//     configurations_->resize( _postoconfs.positions.size(), nbconfs );
//     t_Confs :: const_iterator i_confs = confs.begin();
//     t_Confs :: const_iterator i_confs_end = confs.end();
//     for(size_t j(0); i_confs != i_confs_end; ++i_confs )
//     {
//       PosToConfs :: t_Configurations :: const_iterator i_conf = i_confs->begin();
//       PosToConfs :: t_Configurations :: const_iterator i_conf_end = i_confs->end();
//       for(; i_conf != i_conf_end; ++i_conf, ++j )
//       {
//         __ASSERT( j == nbconfs, "Inconsistent sizes" );
//         for( size_t i(0); i < _postoconfs.positions.size(); ++i )
//           (*configurations_)(i,j) = i_conf->first[i] ? 1: 0;
//       }
//     }
 

      // initializes mapping.
      mapping().init( _strs, confs );
    }

  INCOLLAPSE( opt::ErrorTuple ) :: evaluate() const
  {
    __DEBUGTRYBEGIN
    opt::ErrorTuple error;
    for(size_t n(0); n < mapping().size(); ++n )
      if( not mapping().do_skip(n) ) 
        error += opt::ErrorTuple( mapping().target( n ) - evaluate(n),
                                  mapping().weight( n ) );
    return error;
    __DEBUGTRYEND(, "Error in Collapse::evaluate().\n" )
  }

  INCOLLAPSE( typename COLHEAD::t_Matrix::value_type ) :: evaluate(size_t _n) const 
  {
    __DEBUGTRYBEGIN
    __ASSERT( _n >= mapping().size(), "Inconsistent sizes.\n" )
    namespace bblas = boost::numeric::ublas;
    bblas::range range( mapping().range(_n) );
    types::t_real intermed(0);
    for( bblas::range::const_iterator j( range.begin() ); j != range.end(); ++j )
    {
      const bblas::matrix_column<const t_Configurations> 
        config( *configurations_, *j );
      intermed +=   separables()( config )
                  * mapping().eweight(_n,*j - range.start() );
    }
    return intermed;
    __DEBUGTRYEND(, "Error in Collapse::evaluate().\n" )
  }

  INCOLLAPSE( void ) :: init( t_Separables& _sep )
  {
    separables_ = &_sep; 
    update_.init( _sep );
    regularization().init( _sep );
  }

# undef COLHEAD
# undef INCOLLAPSE
# undef INCOLLAPSE2
} // end of CE namespace.
