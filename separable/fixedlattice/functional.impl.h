//
//  Version: $Id$
//
#include<boost/numeric/ublas/vector_proxy.hpp>
#include<boost/numeric/ublas/matrix_proxy.hpp>
#include<boost/numeric/ublas/operation.hpp>
#include<numeric>

namespace CE
{
  template<class T_TRAITS >
    void Separables<T_TRAITS > :: set_rank_n_size( size_t _rank, size_t _size )
    {
      coefficients_.resize( _rank * t_Mapping::D, _size );
      norms.resize( _rank  ); 
    }
  template<class T_TRAITS >
    size_t Separables<T_TRAITS> :: ranks() const
    {
      __ASSERT( coefficients().size1() % t_Mapping::D, 
                "Inconsistent sizes.\n" )
      return coefficients().size1() / t_Mapping::D; 
    }
  template<class T_TRAITS > template< class T_VECTOR > 
    types::t_real Separables<T_TRAITS> :: operator()( const T_VECTOR &_conf ) const
    {
      namespace bblas = boost::numeric::ublas;
      namespace bl = boost::lambda;
      t_Vector intermed( ranks(), 1e0 );
      t_Policy :: rank_vector( coefficients(), _conf, intermed, bl::_1 *= bl::_2 );
      return bblas::inner_prod( intermed, norms );
    }

  template<class T_TRAITS>
  std::ostream& operator<<( std::ostream& _stream,
                            const Separables<T_TRAITS> &_sep )
  {
    _stream << " Separable Function:";
    typedef Separables<T_TRAITS> t_Separables;
    typedef typename t_Separables :: t_Matrix t_Matrix;
    typename t_Matrix :: const_iterator1 i_row = _sep.coefficients().begin1();
    typename t_Matrix :: const_iterator1 i_row_end = _sep.coefficients().end1();
    for( size_t r(0); i_row != i_row_end; i_row += t_Separables::t_Mapping:: D, ++r )
    {
      _stream << "\n   Rank " << r << ": " << _sep.norms[r] << "\n     ";
      typename t_Matrix :: const_iterator2 i_column = i_row.begin();
      typename t_Matrix :: const_iterator2 i_column_end = i_row.end();
      for( size_t d(0); i_column != i_column_end; ++i_column, ++d )
      {
        _stream << "(";
        typename t_Matrix :: const_iterator1 i_coef = i_column.begin();
        for( size_t i(0); i < t_Separables :: t_Mapping :: D; ++i, ++i_coef )
          _stream  << *i_coef << " ";
        _stream << ") ";
        if( d % 5 == 0 and d ) _stream << "\n     "; 
      }
    }
    return _stream;
  }
  namespace Policy 
  {
#   if defined(POLICYDEF) || defined(POLICYDEF2) || defined(POLICYDEF3) || defined(POLICYDEF4) 
#     error "POLICYDEF macros already exists"
#   endif
#   define POLICYDEF(code) \
      template<class T_MAPPING> code void DimensionMatrix<T_MAPPING> ::
#   define POLICYDEF2(code1, code2) \
      template<class T_MAPPING> code1, code2 void DimensionMatrix<T_MAPPING> ::
#   define POLICYDEF3(code1, code2, code3) \
      template<class T_MAPPING> code1, code2, code3 void DimensionMatrix<T_MAPPING> ::
#   define POLICYDEF4(code1, code2, code3, code4) \
      template<class T_MAPPING> code1, code2, code3, code4 void DimensionMatrix<T_MAPPING> ::
    
    POLICYDEF4( template< class T_COEFS, class T_VECIN, class T_VECOUT, class T_OP > )
      rank_vector( const T_COEFS &_coefs, const T_VECIN &_vecin,
                   T_VECOUT &_vecout, T_OP _op )
      {
        __ASSERT( _coefs.size1() != _vecout.size() * t_Mapping :: D,
                     "Inconsistent sizes: " << _coefs.size1() << " != " 
                  << _vecout.size() << " * " << t_Mapping::D << ".\n" )
        __ASSERT( _vecin.size() != _coefs.size2(), "Inconsistent sizes.\n" )
        typename T_COEFS :: const_iterator2 i_column = _coefs.begin2();
        typename T_COEFS :: const_iterator2 i_column_end = _coefs.end2();
        typename T_VECIN :: const_iterator i_in = _vecin.begin();
        for(; i_column != i_column_end; ++i_column, ++i_in )
        {
          typename T_COEFS :: const_iterator1 i_row = i_column.begin();
          typename T_COEFS :: const_iterator1 i_row_end = i_column.end();
          typename T_VECOUT :: iterator i_out = _vecout.begin();
          for(; i_row != i_row_end; i_row += t_Mapping :: D, ++i_out )
            t_Mapping :: apply( _op, *i_in, i_row, *i_out );
        }
      }

    POLICYDEF2( template< class T_COEFS, class T_NORMS > )
      normalize( T_COEFS &_coefs, T_NORMS &_norms )
      {
        typename T_COEFS :: iterator2 i_column = _coefs.begin2();
        typename T_COEFS :: iterator2 i_column_end = _coefs.end2();
        for(; i_column != i_column_end; ++i_column )
        {
          typename T_COEFS :: iterator1 i_row = i_column.begin();
          typename T_COEFS :: iterator1 i_row_end = i_column.end();
          typename T_NORMS :: iterator i_norm = _norms.begin();
          for(; i_row != i_row_end; i_row += t_Mapping :: D, ++i_norm )
            t_Mapping::normalize( i_row, *i_norm );
        }
      }
    POLICYDEF( template< class T_COEFS > )
      randomize( T_COEFS &_coefs, typename T_COEFS :: value_type _howrandom)
      {
        typename T_COEFS :: iterator2 i_column = _coefs.begin2();
        typename T_COEFS :: iterator2 i_column_end = _coefs.end2();
        for(; i_column != i_column_end; ++i_column )
        {
          typename T_COEFS :: iterator1 i_row = i_column.begin();
          typename T_COEFS :: iterator1 i_row_end = i_column.end();
          for(; i_row != i_row_end; i_row += t_Mapping :: D ) 
            t_Mapping::randomize( i_row, _howrandom );
        }
      }

    POLICYDEF4( template< class T_COEFS, class T_VECIN, class T_OUT, class T_OP > )
      apply_throughout( const T_COEFS &_coefs, const T_VECIN &_vecin, 
                        T_OUT &_out, T_OP _op )
      {
        __ASSERT( _coefs.size1() % t_Mapping :: D != 0, "Inconsistent sizes.\n" )
        std::vector< T_OUT > result( _coefs.size1() / t_Mapping :: D, T_OUT(1) );
        rank_vector( _coefs, _vecin, result,  _op );
        _out = std::accumulate( result.begin(), result.end(), _out );
      }

    POLICYDEF4( template< class T_COEFS, class T_VECIN, class T_VECOUT, class T_OP > )
      apply_to_dim( const T_COEFS &_coefs, const T_VECIN &_vecin,
                    T_VECOUT &_vecout, size_t _d, T_OP _op )
      {
        namespace bblas = boost::numeric::ublas;
        __ASSERT( _coefs.size1() != _vecout.size() * t_Mapping :: D, "Inconsistent sizes.\n" )
        __ASSERT( _vecin.size() != _coefs.size2(), "Inconsistent sizes.\n" )
        __ASSERT( _d >= _coefs.size2(), "Inconsistent input dimension.\n" )
        typedef bblas::matrix_column< T_COEFS > t_Column;
        t_Column column( _coefs, _d );
        typename T_VECIN :: value_type in = _vecin[_d];
        typename t_Column :: const_iterator i_row = column.begin();
        typename t_Column :: const_iterator i_row_end = column.end();
        typename T_VECOUT :: iterator i_out = _vecout.begin();
          for(; i_row != i_row_end; i_row += t_Mapping :: D, ++i_out )
            t_Mapping :: apply( _op, in, i_row, *i_out );
      }
    POLICYDEF4( template< class T_COEFS, class T_VECIN, class T_OUT, class T_OP > )
      apply_to_rank( const T_COEFS &_coefs, const T_VECIN &_vecin,
                     T_OUT &_out, size_t _r, T_OP _op )
      {
        __ASSERT( _r * t_Mapping::D >= _coefs.size1(), "Inconsistent sizes.\n" )
        __ASSERT( _vecin.size() != _coefs.size2(), "Inconsistent sizes.\n" )
        const typename T_COEFS :: const_iterator1 i_row = _coefs.begin1() + _r * t_Mapping::D;
        typename T_COEFS :: const_iterator2 i_column = i_row.begin();
        typename T_COEFS :: const_iterator2 i_column_end = i_row.end();
        typename T_VECIN :: const_iterator i_in = _vecin.begin();
        for(; i_column != i_column_end; ++i_column, ++i_in )
          t_Mapping::apply( _op, *i_in, i_column.begin(), _out );
      }
    POLICYDEF4( template< class T_COEFS, class T_VECIN, class T_OUT, class T_OP > )
      apply_to_dim_n_rank( const T_COEFS &_coefs, const T_VECIN &_vecin,
                           T_OUT &_out, size_t _d, size_t _r, T_OP _op )
      {
        __ASSERT( _vecin.size() != _coefs.size2(), "Inconsistent sizes.\n" )
        __ASSERT( _d >= _coefs.size2(), "Inconsistent input dimension.\n" )
        typename T_VECIN :: value_type in = _vecin[_d];
        t_Mapping :: apply( _op, in, 
                            (_coefs.begin2() + _d).begin() + _r * t_Mapping::D,
                            _out );
      }
#   undef POLICYDEF
#   undef POLICYDEF2
#   undef POLICYDEF3
#   undef POLICYDEF4
  } // end of Policy namespace
} // end of CE namespace.
