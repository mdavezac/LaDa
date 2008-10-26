//
//  Version: $Id$
//
#include<boost/numeric/ublas/vector_proxy.hpp>
#include<boost/numeric/ublas/matrix_proxy.hpp>
#include<boost/numeric/ublas/operation.hpp>
#include<boost/lambda/bind.hpp>

namespace CE
{
  namespace Policy 
  {
    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      typename T_SEPARABLES :: t_Vector :: value_type
        LowMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS>
          :: factor( size_t _kv, size_t _r, size_t _d,
                     const t_Separables &_sep,
                     const t_Configurations &_confs ) const
          {
            namespace bl = boost::lambda;
            namespace bblas = boost::numeric::ublas;
         
            typename t_Vector :: value_type result;
            bblas::matrix_column< const t_Configurations > config( _confs, _kv );
            t_Separables :: t_Policy
                         :: apply_to_dim_n_rank( _sep.coefficients(),
                                                 config, result, _d, _r,
                                                 bl::_1 = bl::_2 );
            return not Fuzzy::is_zero( result ) ?
                     scales_(_r, _kv) / result:
                     factor_from_scratch( _kv, _r, _d, _sep, _confs );
          }
    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      typename T_SEPARABLES :: t_Vector :: value_type
        LowMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS>
          :: factor_from_scratch( size_t _kv, size_t _r, size_t _d,
                                  const t_Separables &_sep,
                                  const t_Configurations &_confs ) const
          {
            namespace bl = boost::lambda;
            namespace bblas = boost::numeric::ublas;
         
            typename t_Vector::value_type result(1);
            size_t d(0);
            bblas::matrix_column< const t_Configurations > config( _confs, _kv );
            t_Separables :: t_Policy :: apply_to_rank
            ( 
              _sep.coefficients(), config, result, _r,
              (
                bl::if_then( bl::var(d) != bl::constant(_d), bl::_1 *= bl::_2  ),
                ++bl::var(d)
              )
            );
            return result;
          }

    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      void LowMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS>
        :: operator()( const t_Separables &_sep,
                       const t_Configurations &_confs )
        {
          namespace bl = boost::lambda;
          namespace bblas = boost::numeric::ublas;
          if( not (     scales_.size2() == _confs.size2() 
                    and scales_.size1() == _sep.ranks() )  )
              scales_.resize( _sep.ranks(), _confs.size2() );
          for( size_t i(0); i < _confs.size2(); ++i )
          {
            bblas::matrix_column<const t_Configurations> conf( _confs, i );
            bblas::matrix_column<t_CMatrix> scaling( scales_, i );
            std::fill( scaling.begin(), scaling.end(), typename t_Matrix::value_type(1) );
            t_Separables::t_Policy::rank_vector
            ( 
              _sep.coefficients(), conf, scaling,
              bl::_1 *= bl::_2
            );
          }
        }

    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      void HighMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS> 
        :: operator()( const t_Separables &_sep,
                       const t_Configurations &_confs )
        {
          // First updates stuff split along dimensions.
          namespace bl = boost::lambda;
          namespace bblas = boost::numeric::ublas;

          if( not (     scales_.size2() == _confs.size2() 
                    and scales_.size1() == _sep.ranks() )  )
          {
            scales_.resize( _sep.ranks(), _confs.size2() );
            dimsplit_.resize( _confs.size2() ); 
            typename std::vector< t_CMatrix > :: iterator i_split = dimsplit_.begin();
            typename std::vector< t_CMatrix > :: iterator i_split_end = dimsplit_.end();
            for(; i_split != i_split_end; ++i_split )
              i_split->resize( _sep.ranks(), _sep.dimensions() );
          }

          // First updates stuff split along dimensions.
          typename std::vector< t_CMatrix > :: iterator i_split = dimsplit_.begin();
          typename std::vector< t_CMatrix > :: iterator i_split_end = dimsplit_.end();
          for( size_t i(0); i_split != i_split_end; ++i, ++i_split )
          {
            bblas::matrix_column< const t_Configurations > config( _confs, i );
            for( size_t d(0); d < _sep.dimensions(); ++d ) 
              for(size_t r(0); r < _sep.ranks(); ++r )
                t_Separables :: t_Policy
                             :: apply_to_dim_n_rank( _sep.coefficients(),
                                                     config, (*i_split)(r,d), d, r,
                                                     bl::_1 = bl::_2 );
          }
       
          // Then  updates scales_.
          update_scales( _sep );
        }
    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      void HighMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS>
        :: update_scales( const t_Separables &_sep )
        {
          namespace bl = boost::lambda;
          namespace bblas = boost::numeric::ublas;
          typedef typename std::vector< t_CMatrix > :: const_iterator t_cit;
          t_cit i_split = dimsplit_.begin();
          t_cit i_split_end = dimsplit_.end();
          for( size_t i(0); i_split != i_split_end; ++i, ++i_split )
          {
            for(size_t r(0); r < _sep.ranks(); ++r )
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
      void HighMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS> 
        :: operator()( size_t _dim,
                       const t_Separables &_sep,
                       const t_Configurations &_confs )
        {
          // First updates stuff split along dimensions.
          namespace bl = boost::lambda;
          namespace bblas = boost::numeric::ublas;
          // First updates stuff split along dimensions.
          typename std::vector< t_CMatrix > :: iterator i_split = dimsplit_.begin();
          typename std::vector< t_CMatrix > :: iterator i_split_end = dimsplit_.end();
          for( size_t i(0); i_split != i_split_end; ++i, ++i_split )
          {
            bblas::matrix_column< const t_Configurations > config( _confs, i );
            for(size_t r(0); r < _sep.ranks(); ++r )
              t_Separables :: t_Policy
                           :: apply_to_dim_n_rank( _sep.coefficients(),
                                                   config, (*i_split)(r, _dim),
                                                   _dim, r, bl::_1 = bl::_2 );
          }
 
          // Then  updates scales_.
          update_scales( _sep );
        }
    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      typename T_SEPARABLES :: t_Vector :: value_type
        HighMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS>
          :: factor( size_t _kv, size_t _r, size_t _d,
                     const t_Separables &_sep,
                     const t_Configurations &_confs ) const
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
                     row.begin(), row.end(), 1e0,
                     bl::ret< typename t_CMatrix :: value_type >
                     ( 
                       bl::if_then_else_return
                       ( 
                         bl::var(d) != bl::constant(_d),
                         bl::_1 * bl::_2,
                         bl::_1
                       )
                     )
                   );
          }
        
    template< class T_SEPARABLES > template< class T_MATRIX, class T_VECTOR >
      void BadRegularization<T_SEPARABLES> 
        :: operator()( const t_Separables& _sep, T_MATRIX &_A,
                       T_VECTOR &_b, size_t _dim ) const
        {
          if( Fuzzy::is_zero( lambda ) ) return; 
          typedef typename t_Separables :: t_Vector :: const_iterator t_cit;
          t_cit i_norm = _sep.norms.begin();
          t_cit i_norm_end = _sep.norms.end();
          for( size_t j(0); i_norm != i_norm_end; ++i_norm )
          {
            typedef typename T_MATRIX::value_type t_Type;
            t_Type factor( lambda  * (*i_norm) * (*i_norm) );
            for( size_t i(0); i < t_Separables :: t_Mapping :: D; ++i, ++j )
              _A( j, j ) +=  factor * t_Type( t_Separables::t_Mapping::norm( i ) );
          }
        }
    template< class T_SEPARABLES > template< class T_MATRIX, class T_VECTOR >
      void Regularization<T_SEPARABLES> 
        :: operator()( const t_Separables &_sep, T_MATRIX &_A,
                       T_VECTOR &_b, size_t _dim ) const
        {
          if( Fuzzy::is_zero( lambda ) ) return; 
          const size_t D( t_Separables :: t_Mapping :: D );
          typedef typename t_Separables :: t_Vector :: const_iterator t_cit;
          t_cit i_norm = _sep.norms.begin();
          t_cit i_norm_end = _sep.norms.end();
          for( size_t j(0); i_norm != i_norm_end; ++i_norm, j += D )
          {
            typedef typename T_MATRIX::value_type t_Type;
            t_Type factor( lambda ); //* (*i_norm) * (*i_norm) );
            for( size_t i(0); i <  D; ++i )
              for( size_t k(0); k < D; ++k )
                if( i == k ) _A( j + i, j + k ) += factor *( 1.e0 - 1.e0/t_Type(D) );
                else  _A( j + i, j + k ) += factor *( -1.e0/t_Type(D) );
          }
        }

//   template< class T_TRAITS > 
//     void Initialization :: init< 


 } // end of Policy namespace
} // end of CE namespace.
