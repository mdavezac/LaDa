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
          :: factor( size_t _kv, size_t _r, size_t _d ) const
          {
            namespace bl = boost::lambda;
            namespace bblas = boost::numeric::ublas;
         
            typename t_Vector :: value_type result;
            bblas::matrix_column< const t_Configurations >
              config( *configurations_, _kv );
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
            bblas::matrix_column< const t_Configurations >
              config( *configurations_, _kv );
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
        __ASSERT( scales_.size2() != configurations_->size2(),
                  "Inconsistent sizes.\n" )
        for( size_t i(0); i < configurations_->size2(); ++i )
        {
          bblas::matrix_column<const t_Configurations> conf( *configurations_, i );
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
          scales_.resize( separables_->ranks(), configurations_->size2() ); 
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
          bblas::matrix_column< const t_Configurations > 
            config( *configurations_, i );
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
        __ASSERT( scales_.size2() != configurations_->size2(),
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
          bblas::matrix_column< const t_Configurations >
            config( *configurations_, i );
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
          dimsplit_.resize( configurations_->size2() ); 
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
} // end of CE namespace.
