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
      __ASSERT( dim >= configurations_.size1(), "Inconsistent sizes.\n" )
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
    
        _A += mapping().weight(i) * bblas::outer_prod( X, X ); 
        _b += mapping().weight(i) * mapping().target(i) * X;
      }
    }

  INCOLLAPSE(void) :: create_X( size_t _i, t_Vector &_out )
  {
    namespace bblas = boost::numeric::ublas;
    __ASSERT( dof() != _out.size(), "Incompatible sizes.\n" )
    // Create matrix range including only equivalent configurations.
    typedef const bblas::matrix_range< t_iMatrix > t_Range;
    const bblas::range equivrange( mapping().range( _i ) );
    const bblas::range dimrange( dim, dim+1 );
    t_Range equivconfs( configurations_, dimrange, equivrange );

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
        bblas::vector_range< t_Vector > vecr( _out, range );
        t_Separables::t_Mapping::add_tovec( *i_conf, vecr, scalar );
      }
  }

  INCOLLAPSE(void) :: update_all()
  {
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
      mapping().init( _strs, confs );
    }

  INCOLLAPSE( opt::ErrorTuple ) :: evaluate() const
  {
    namespace bblas = boost::numeric::ublas;
    opt::ErrorTuple error;
    for(size_t n(0); n < mapping().size(); ++n )
    {
      if( mapping().do_skip(n) ) continue;
      bblas::range range( mapping().range(n) );
      types::t_real intermed(0);
      for( bblas::range::const_iterator j( range.begin() ); j != range.end(); ++j )
      {
        const bblas::matrix_column<const t_iMatrix> config( configurations_, *j );
        intermed +=   separables()( config )
                    * mapping().eweight(n,*j - range.start() );
      }
      error += opt::ErrorTuple( mapping().target(n) - intermed, mapping().weight(n) );
    }
    return error;
  }

  INCOLLAPSE( void ) :: init( t_Separables& _sep )
  {
    separables_ = &_sep; 
    update_.init( _sep );
    regularization().init( _sep );
  }

  namespace Policy 
  {
    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      typename T_SEPARABLES :: t_Vector :: value_type
        LowMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS>
          :: factor( size_t _kv, size_t _r, size_t _d ) const
          {
            namespace bl = boost::lambda;
            namespace bblas = boost::numeric::ublas;
         
            typename t_Vector :: value_type result;
            bblas::matrix_column< const t_iMatrix > config( configurations_, _kv );
            t_Separables :: t_Policy
                         :: apply_to_dim_n_rank( separables_->coefficients(),
                                                 config, result, _d, _r,
                                                 bl::_1 = bl::_2 );
            return not Fuzzy::is_zero( result ) ?
                     scales_(_r, _kv) / result:
                     factor_from_scratch( _kv, _r, _d );
          }
    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      typename T_SEPARABLES :: t_Vector :: value_type
        LowMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS>
          :: factor_from_scratch( size_t _kv, size_t _r, size_t _d ) const
          {
            namespace bl = boost::lambda;
            namespace bblas = boost::numeric::ublas;
         
            typename t_Vector::value_type result(1);
            size_t d(0);
            bblas::matrix_column< const t_iMatrix > config( configurations_, _kv );
            t_Separables :: t_Policy :: apply_to_rank
            ( 
              separables_->coefficients(), config, result, _r,
              (
                bl::if_then( bl::var(d) != bl::constant(_d), bl::_1 *= bl::_2  ),
                ++bl::var(d)
              )
            );
            return result;
          }

    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      void LowMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS> :: operator()()
      {
        namespace bl = boost::lambda;
        namespace bblas = boost::numeric::ublas;
        __ASSERT( scales_.size2() != configurations_.size2(),
                  "Inconsistent sizes.\n" )
        for( size_t i(0); i < configurations_.size2(); ++i )
        {
          bblas::matrix_column<const t_iMatrix> conf( configurations_, i );
          bblas::matrix_column<t_CMatrix> scaling( scales_, i );
          std::fill( scaling.begin(), scaling.end(), typename t_Matrix::value_type(1) );
          t_Separables::t_Policy::rank_vector
          ( 
            separables_->coefficients(), conf, scaling,
            bl::_1 *= bl::_2
          );
        }
      }
    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      void LowMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS>
        :: init( const t_Separables& _sep )
        {
          separables_ = &_sep;
          scales_.resize( separables_->ranks(), configurations_.size2() ); 
        }

    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      void HighMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS> :: operator()()
      {
        // First updates stuff split along dimensions.
        namespace bl = boost::lambda;
        namespace bblas = boost::numeric::ublas;
        // First updates stuff split along dimensions.
        typename std::vector< t_CMatrix > :: iterator i_split = dimsplit_.begin();
        typename std::vector< t_CMatrix > :: iterator i_split_end = dimsplit_.end();
        for( size_t i(0); i_split != i_split_end; ++i, ++i_split )
        {
          bblas::matrix_column< const t_iMatrix > config( configurations_, i );
          for( size_t d(0); d < separables_->dimensions(); ++d ) 
            for(size_t r(0); r < separables_->ranks(); ++r )
              t_Separables :: t_Policy
                           :: apply_to_dim_n_rank( separables_->coefficients(),
                                                   config, (*i_split)(r,d), d, r,
                                                   bl::_1 = bl::_2 );
        }

        // Then  updates scales_.
        update_scales();
      }
    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      void HighMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS> :: update_scales()
      {
        namespace bl = boost::lambda;
        namespace bblas = boost::numeric::ublas;
        __ASSERT( scales_.size2() != configurations_.size2(),
                  "Inconsistent sizes.\n" )
        typedef typename std::vector< t_CMatrix > :: const_iterator t_cit;
        t_cit i_split = dimsplit_.begin();
        t_cit i_split_end = dimsplit_.end();
        for( size_t i(0); i_split != i_split_end; ++i, ++i_split )
        {
          for(size_t r(0); r < separables_->ranks(); ++r )
          {
            bblas::matrix_row< const t_CMatrix > row( *i_split, r );
            scales_(r,i) = std::accumulate
                           (
                             row.begin(), row.end(),
                             typename t_CMatrix::value_type(1),
                             bl::_1 * bl::_2
                           );
          }
        }
      }
    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      void HighMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS> :: operator()(size_t _dim )
      {
        // First updates stuff split along dimensions.
        namespace bl = boost::lambda;
        namespace bblas = boost::numeric::ublas;
        // First updates stuff split along dimensions.
        typename std::vector< t_CMatrix > :: iterator i_split = dimsplit_.begin();
        typename std::vector< t_CMatrix > :: iterator i_split_end = dimsplit_.end();
        for( size_t i(0); i_split != i_split_end; ++i, ++i_split )
        {
          bblas::matrix_column< const t_iMatrix > config( configurations_, i );
          for(size_t r(0); r < separables_->ranks(); ++r )
            t_Separables :: t_Policy
                         :: apply_to_dim_n_rank( separables_->coefficients(),
                                                 config, (*i_split)(r, _dim),
                                                 _dim, r, bl::_1 = bl::_2 );
        }

        // Then  updates scales_.
        update_scales();
      }
    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      typename T_SEPARABLES :: t_Vector :: value_type
        HighMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS>
          :: factor( size_t _kv, size_t _r, size_t _d ) const
          {
            typename t_Matrix :: value_type result( dimsplit_[_kv](_r,_d) );
            return not Fuzzy::is_zero( result ) ?
                     scales_(_r, _kv) / result:
                     factor_from_scratch( _kv, _r, _d );
          }
    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      typename T_SEPARABLES :: t_Vector :: value_type
        HighMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS>
          :: factor_from_scratch( size_t _kv, size_t _r, size_t _d ) const
          {
            namespace bl = boost::lambda;
            namespace bblas = boost::numeric::ublas;
         
            typename bblas::matrix_row<const t_CMatrix> row( dimsplit_[_kv], _r );
            size_t d(0);
            return std::accumulate
                   (
                     row.begin(), row.end(), d,
                     bl::if_then_else_return
                     ( 
                       bl::var(d) != bl::constant(_d),
                       bl::_1 * bl::_2,
                       bl::_1
                     )
                   );
          }
    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      void HighMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS>
        :: init( const t_Separables& _sep )
        {
          namespace bl = boost::lambda;
          t_Base :: init( _sep );
          dimsplit_.resize( configurations_.size2() ); 
          typename std::vector< t_CMatrix > :: iterator i_split = dimsplit_.begin();
          typename std::vector< t_CMatrix > :: iterator i_split_end = dimsplit_.end();
          for(; i_split != i_split_end; ++i_split )
            i_split->resize( separables_->ranks(), separables_->dimensions() );
        }
    template< class T_SEPARABLES > template< class T_MATRIX, class T_VECTOR >
      void Regularization<T_SEPARABLES> 
        :: operator()( T_MATRIX &_A, T_VECTOR &_b, size_t _dim )
        {
          if( Fuzzy::is_zero( lambda ) ) return; 
          typedef typename t_Separables :: t_Vector :: const_iterator t_cit;
          t_cit i_norm = separables_->norms.begin();
          t_cit i_norm_end = separables_->norms.end();
          for( size_t j(0); i_norm != i_norm_end; ++i_norm )
          {
            typedef typename T_MATRIX::value_type t_Type;
            t_Type factor( lambda  * (*i_norm) * (*i_norm) );
            for( size_t i(0); i < t_Separables :: t_Mapping :: D; ++i, ++j )
              _A( j, j ) +=  factor * t_Type( t_Separables::t_Mapping::norm( i ) );
          }
        }

 } // end of Policy namespace
# undef COLHEAD
# undef INCOLLAPSE
# undef INCOLLAPSE2
} // end of CE namespace.
