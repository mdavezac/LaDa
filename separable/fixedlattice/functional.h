//
//  Version: $Id$
//
#ifndef _CE_SEPARABLE_H_
#define _CE_SEPARABLE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/matrix.hpp>

#include <opt/types.h>
#include <opt/debug.h>

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


  //! \brief A sum of separable functions for a fixed-lattice.
  //! \param T_MAPPING specifies how to map from configuration to coefficients
  //!                  elements. This mapping takes care of representing the
  //!                  number of species adequately. There are at present two
  //!                  types: VectorPlus where each type is mapped to vectors
  //!                  (1,0,..), (0,1,...), ..., and a VectorDiff, for
  //!                  regularization, for which each element is mapped to a
  //!                  constant vector (1,1,1,...) and then subsequently to
  //!                  (0,1,0,..,.) vectors. 
  //! \details Separables::coefficients are arrange in a matrix where each row
  //!          is a different rank of the sum of seps. The colums are organized
  //!          along direction and inner basis for direction. The inner basis
  //!          is the faster running index.
  template< class T_MAPPING >
    class Separables
    {
      public:
        //! Type of mapping used to go from conf to coefs.
        typedef T_MAPPING t_Mapping;
        //! Type of boost matrices.
        typedef boost::numeric::ublas::matrix<types::t_real> t_Matrix;
        //! Type of boost vectors.
        typedef boost::numeric::ublas::vector<types::t_real> t_Vector;

        //! \brief Coefficients of the separable functions. 
        //! \details Rows denote ranks, and columns are indexed according to
        //!          the dimension of the separable function and the family of
        //!          functions of each separable function.
        t_Matrix coefficients;
        //! Holds norms of each rank.
        t_Vector norms;

        //! Constructor.
        Separables() {}
        //! Copy constructor.
        Separables   ( const Separables<t_Mapping> &_c )
                   : coefficients( _c.coefficients ), norms( _c.norms ) {}
        //! Destructor.
        ~Separables() {}

        //! Returns the value of the separable function evaluated at \a _conf.
        template< class T_VECTOR >
        types::t_real operator()( const T_VECTOR &_conf ) const
        {
          namespace bblas = boost::numeric::ublas;
          t_Vector intermed( coefficients.size1() );
          std::fill( intermed.begin(), intermed.end(), 0e0 );
          details::rank_vector< t_Matrix, t_Vector, T_VECTOR, t_Mapping>
                              ( coefficients, _conf, intermed );
          return bblas::inner_prod( intermed, norms );
        }
        //! Sets ranks and sizes.
        void set_rank_n_size( size_t _rank, size_t _size )
         { coefficients.resize( _rank, _size * t_Mapping :: D ); } 
    };

  template<class T_MAPPING >
  std::ostream& operator<<( std::ostream& _stream, const Separables<T_MAPPING> &_sep )
  {
    _stream << " Separable Function:\n";
    for( size_t i(0); i < _sep.coefficients.size1(); ++i )
    {
      _stream << "   Rank " << i << ": " << _sep.norms[i] << "\n     `";
      for( size_t j(0); j < _sep.coefficients.size2(); ++j )
      {
        _stream << "(";
        for( size_t d(0); d >= T_MAPPING::D; ++d )
         _stream  << _sep.coefficients( i, j * T_MAPPING::D + d ) << " ";
        _stream << ") ";
        if( j % 5 == 0 and j ) _stream << "\n     "; 
      }
      if( _sep.coefficients.size2() % 5 != 0 ) _stream << "\n";
    }
    return _stream;
  }
}

#include "functional.impl.h"

#endif
