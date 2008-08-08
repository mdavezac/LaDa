//
//  Version: $Id$
//
#ifndef _SEPARABLE_CE_H_
#define _SEPARABLE_CE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <vector>

#include <opt/types.h>

#include "separable.h"
#include "boolean.h"

namespace CE
{

  // Forward Declarations.
  //! \cond
  namespace details 
  {
    template< class T_MATRIX, class T_VECTOR1, class T_VECTOR2, class T_MAPPING > 
      void rank_vector( const T_MATRIX &_mat, const T_VECTOR1 &_in,
                        T_VECTOR2 &_out, const T_MAPPING &_mapping );
 
    template< class T_VECTOR1, class T_VECTOR2, class T_MAPPING > 
      typename T_VECTOR1::value_type conf_coef( const T_VECTOR1 &_coef,
                                                const T_VECTOR2 &_conf,
                                                const T_MAPPING &_mapping );

    template< class T_MATRIX1, class T_MATRIX2, class T_VECTOR, class T_MAPPING > 
      void allconfs_rank_vector( const T_MATRIX1 &_mat, const T_MATRIX2 &_in,
                                 T_VECTOR &_out, const T_MAPPING &_mapping );

    template< class T_MATIN, class T_MATOUT, class T_VEC, class T_MAPPING > 
      void rank_dim_matrix( const T_MATIN &_in, const T_VEC &_vec,
                            T_MATOUT &_out, const T_MAPPING &_mapping );

    template< class T_VECTOR1, class T_VECTOR2, class T_VECTOR3, class T_MAPPING > 
      typename T_VECTOR1::value_type conf_coef_vector( const T_VECTOR1 &_coef,
                                                       const T_VECTOR2 &_conf,
                                                       const T_VECTOR3 &_out,
                                                       const T_MAPPING &_mapping );

    template< size_t > class VectorPlus;
    template< size_t > class VectorDiff;
  }


  template< class T_MAPPING, class T_NORMALIZATION >
    class Separables
    {
      public:
        //! Type of mapping used to go from conf to coefs.
        typedef T_MAPPING t_Mapping;
        //! Type of normalization used.
        typedef T_NORMALIZATION t_Normalization;
        //! Type of boost matrices.
        boost::numeric::ublas::matrix<types::t_real> t_Matrix;
        //! Type of boost vectors.
        boost::numeric::ublas::matrix<types::t_real> t_Vector;

        //! Returns the value of the separable function evaluated at \a _conf.
        template< class T_VECTOR >
        types::t_real operator()( const T_VECTOR &_conf ) const
        {
          t_Vector intermed( coefficients.size1() );
          std::fill( intermed.begin(), intermed.end(), 0e0 );
          details::rank_vector< t_Matrix, t_Vector, T_VECTOR, t_Mapping>
                              ( coefficients, _conf, intermed );
          return bblas::inner_prod( intermed, norms );
        }

      public:
        //! \brief Coefficients of the separable functions. 
        //! \details Rows denote ranks, and columns are indexed according to
        //!          the dimension of the separable function and the family of
        //!          functions of each separable function.
        t_Matrix coefficients;
        //! Holds norms of each rank.
        t_Vector norms;
    };

    //! \brief Returns the sum of the separable function \a _seps evaluated
    //!        over each column of \a _confs.
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


  namespace details 
  {
    template< class T_MATRIX, class T_VECTOR, size_t D, class T_MAPPING > 
      void rank_vector( const T_MATRIX &_mat, const T_VECTOR &_in,
                        T_VECTOR &_out, const T_MAPPING &_mapping )
      {
        namespace bblas = boost::numeric::ublas;
        __ASSERT( _in.size() * D != mat.size2(),
                  "Input vector and matrix of incompatible size.\n" )
        __ASSERT( _out.size() != mat.size1(),
                  "Output vector and matrix of incompatible size.\n" )
        t_Vector :: iterator i_out = _out.begin();
        t_Vector :: iterator i_out_end = _out.end();
        for(size_t i(0); i_out != i_out_end; ++i_out, ++i )
        {
          bblas::matrix_row< T_MATRIX > row( _mat(i) );
          *i_out += conf_coef<T_VECTOR, T_MAPPIND, D>( row, _in, _mapping );
        }
      }
 
    template< class T_VECTOR, class T_MAPPING, size_t D > 
      typename T_VECTOR::value_type conf_coef( const T_VECTOR &_coef,
                                               T_VECTOR &_conf,
                                               const T_MAPPING &_mapping )
      {
        namespace bblas = boost::numeric::ublas;
        __ASSERT( _coef.size() * D != conf.size(),
                  "Coef vector and Conf vector of incompatible size.\n" )
        typename T_VECTOR :: value_type result(0);
        t_Vector :: iterator i_conf = _conf.begin();
        t_Vector :: iterator i_conf_end = _conf.end();
        t_Vector :: iterator i_coef = _coef.end();
        for(; i_conf != i_conf_end; ++i_conf, i_coef += D )
          result += *( i_coef + mapping(  *i_conf ) );
      }
  }
}
