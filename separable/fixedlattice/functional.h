//
//  Version: $Id$
//
#ifndef _SEPARABLE_CE_H_
#define _SEPARABLE_CE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/blas/numeric/vector.hpp>
#include<boost/blas/numeric/matrix.hpp>

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
        boost::numeric::ublas::matrix<types::t_real> t_Matrix;
        //! Type of boost vectors.
        boost::numeric::ublas::matrix<types::t_real> t_Vector;

        //! Constructor.
        Separables() {}
        //! Copy constructor.
        Separables   ( const Separables<t_Mapping, t_Normalization> &_c )
                   : coefficients( _c.coefficients ), norm( _c.norm ) {}
        //! Destructor.
        ~Separables() {}

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

  //! Collapse functor for fitting CE::Separables
  template< class T_SEPARABLES, class T_MAPPING,
            class T_NORMALIZATION = typename T_SEPARABLES :: t_Mapping >
    class Collapse
    {
      public:
        //! Type of the separable function.
        typedef T_SEPARABLES t_Separables;
        //! Type of the mapping function from structures to targets.
        typedef T_MAPPING t_Mapping;
        //! Type of normalization used.
        typedef T_NORMALIZATION t_Normalization;
        //! Type of the Normalization function.
        typedef typename t_Separable :: t_Normalization t_Normalization;
        //! Type of the matrices.
        typedef typename t_Separable :: t_Matrix t_Matrix;
        //! Type of the vectors.
        typedef typename t_Separable :: t_Vector t_Vector;

        //! \brief The mapping from target values to symetrically equivalent
        //!        structures.
        t_Mapping mapping;

        //! Constructor.
        Collapse() : dim(0) {}
        //! Destructor.
        ~Collapse() {}

        //! Creates the fitting matrix and target vector.
        void operator()( t_Matrix &_a, t_Vector &_b,
                         const t_Vector&, const t_Matrix &_coef,
                         types::t_unsigned &_dim );
          { dim = _dim; create_A_n_b( _A, _b, _coef ); }

        //! Updates the scales vector and  normalizes.
        void update_all( t_Matrix &_coefs )
        //! Updates the scales vector, should update only one dim.
        void update( types::t_unsigned _d, t_Matrix &_coefs )
          { update_all( _coefs ); }
        //! Does nothing.
        void reset() {}

        //! Initializes collapse functor.
        template< class T_STRUCTURES >
        void init( const T_STRUCTURES& _strs, const std::string& _bdesc );

      protected:
        //! Creates the _A and _b matrices for fitting.
        void create_A_n_b( t_Matrix &_A, t_Matrix &_b, const t_Vector &_scal );
        //! Finds scaling factor for that conf, collapsed dimension, and rank.
        typename t_Vector::value_type factor( t_Matrix &_coef, size_t _kv, size_t _r );
        //! \brief Creates an X vector for fitting, for a single rank.
        //! \details Before the scaling is done, all ranks equal.
        void create_X( size_t _i, t_Vector &_out, const t_Matrix &_coefs  );


        //! The configurations, arranged in columns.
        t_Matrix configurations;
        //! Holds current dimension being fitted.
        size_t dim;
        //! holds the sep. func. split up by rank and confs.
        t_Matrix scales;
        //! Holds the normalizations. 
        t_Vector norm_vec;
    };
}
