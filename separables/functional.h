#ifndef _CE_SEPARABLE_H_
#define _CE_SEPARABLE_H_

#include "LaDaConfig.h"

#include<boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/matrix.hpp>

#include <misc/types.h>
#include <opt/debug.h>
#include <opt/indirection.h>

namespace LaDa
{
  namespace CE
  {

    // Forward Declarations.
    //! \cond
    namespace Mapping 
    {
      template< size_t > class VectorDiff;
    }
    namespace Policy 
    {
      template< class > class DimensionMatrix;
      class DirectCoefficients;
      class MatrixRangeCoefficients;
    }
    //! \endcond
  }

  namespace Traits
  {
    namespace CE
    {
      template< class T_MAPPING = ::LaDa::CE::Mapping::VectorDiff<2>, 
                class T_POLICY = ::LaDa::CE::Policy::DimensionMatrix< T_MAPPING >,
                class T_COEFFICIENTS = ::LaDa::CE::Policy::DirectCoefficients, 
                class T_VECTOR = boost::numeric::ublas::vector
                        < typename T_COEFFICIENTS :: t_Matrix :: value_type > >
       struct Separables
       {
         //! Type of mapping.
         typedef T_MAPPING t_Mapping;
         //! Type of policy.
         typedef T_POLICY t_Policy;
         //! Type of the coefficient interface.
         typedef T_COEFFICIENTS t_Coefficients;
         //! Type of matrix.
         typedef typename t_Coefficients :: t_Matrix t_Matrix;
         //! Type of vector.
         typedef T_VECTOR t_Vector;
         //! Allows rebinding of the separables function traits.
         template< class TT_MAPPING, class TT_POLICY, 
                   class TT_COEFFICIENTS, class TT_VECTOR >
          struct rebind
          {
            //! new type.
            typedef Separables< TT_MAPPING, TT_POLICY,
                                TT_COEFFICIENTS, TT_VECTOR > type;
          };
       };

      template< class T_SEPARABLES >  struct SeparablesWithMatrixRange;
    } // end of CE namespace 
  } // end of Traits namespace

  namespace CE
  {
    //! \brief A sum of separable functions for a fixed-lattice.
    //! \param T_MAPPING specifies how to map from configuration to coefficients
    //!                  elements. This mapping takes care of representing the
    //!                  number of species adequately. There are at present two
    //!                  types: VectorPlus where each type is mapped to vectors
    //!                  (1,0,..), (0,1,...), ..., and a VectorDiff, for
    //!                  regularization, for which each element is mapped to a
    //!                  constant vector (1,1,1,...) and then subsequently to
    //!                  (0,1,0,..,.) vectors. 
    //! \details Separables::coefficients_ are arrange in a matrix where each column 
    //!          is a different dimension. The row are organized according to
    //!          rank and inner basis (for that column's dimension) with the
    //!          latter the fastest runnning index. Separables::coefficients_
    //!          provides an indirection to this matrix, such that it can be a
    //!          matrix range, a pointer to a matrix, or a direct instantiation
    //!          of a matrix, depending on the type of T_TRAITS::t_Coefficients.
    template< class T_TRAITS >
      class Separables
      {
        public:
          //! bullshit
          const static bool has_norms = true;
          //! Allows rebinding of the separables function.
          template< class TT_TRAITS > struct rebind
          {
            //! new separable type.
            typedef Separables< TT_TRAITS > type;
          };
          //! Type of traits for this class.
          typedef T_TRAITS t_Traits;
          //! Type of mapping used to go from conf to coefs.
          typedef typename t_Traits :: t_Mapping t_Mapping;
          //! Type of boost matrices.
          typedef typename t_Traits :: t_Matrix t_Matrix;
          //! Type of boost vectors.
          typedef typename t_Traits :: t_Vector t_Vector;
          //! Type of the policy class.
          typedef typename t_Traits :: t_Policy t_Policy;
          //! Type of the coefficient interface.
          typedef typename t_Traits :: t_Coefficients t_Coefficients;

          //! Holds norms of each rank.
          t_Vector norms;

          //! Constructor.
          Separables() {}
          //! Copy constructor.
          Separables   ( const Separables &_c )
                     : coefficients_( _c.coefficients_ ), norms( _c.norms ) {}
          //! Destructor.
          ~Separables() {}

          //! Returns the value of the separable function evaluated at \a _conf.
          template< class T_VECTOR >
            types::t_real operator()( const T_VECTOR &_conf ) const;
          //! Sets ranks and sizes.
          void set_rank_n_size( size_t _rank, size_t _size );
          //! Normalizes coefficients.
          void normalize() { t_Policy::normalize( coefficients(), norms ); }
          //! Randomizes coefficients.
          void randomize( typename t_Vector :: value_type _howrandom )
          {
            std::fill( norms.begin(), norms.end(), typename t_Vector::value_type(1) );
            t_Policy::randomize( coefficients(), _howrandom ); 
            normalize();
          }

          //! Returns number of ranks.
          size_t ranks() const;
          //! Returns number of dimensions;
          size_t dimensions() const { return coefficients().size2(); }
          //! Returns the number of degrees of liberty (per dimension).
          size_t dof() const { return coefficients().size1(); }

          //! Returns a reference to the coefficients.
          t_Matrix& coefficients() { return coefficients_(); }
          //! Returns a constant reference to the coefficients.
          const t_Matrix& coefficients() const { return coefficients_(); }

          //! Allows manipulation of the coefficients' interface itself.
          t_Coefficients& coefficients_interface() { return coefficients_; }
          //! Allows manipulation of the coefficients' interface itself.
          const t_Coefficients& coefficients_interface() const { return coefficients_; }


          //! Serializes a separable function.
          template<class ARCHIVE> void serialize( ARCHIVE & _ar, 
                                                  const unsigned int _version);

        protected:
          //! \brief Coefficients of the separable functions. 
          //! \details Rows denote ranks, and columns are indexed according to
          //!          the dimension of the separable function and the family of
          //!          functions of each separable function.
          t_Coefficients coefficients_;
      };

    //! Prints out the separable to a stream.
    template<class T_TRAITS >
    std::ostream& operator<<( std::ostream& _stream,
                              const Separables<T_TRAITS> &_sep );

    //! Holds policies for fixed lattice separables.
    namespace Policy
    {
      //! \brief Holds operations for dimensional matrices.
      //! \details Dimensional matrices in fixed-lattice separable functions are
      //!          organized such that a complete dimension can be recovered as a
      //!          vector from each column of the coefficient matrix.
      template<class T_MAPPING> 
      struct DimensionMatrix
      {
        //! Helps retype the policy.
        template< class T_MAP > struct rebind
        {
          //! Type of the rebound policy.
          typedef DimensionMatrix<T_MAP> type; 
        };
        //! type of the mapping.
        typedef T_MAPPING t_Mapping;
        //! \brief Returns a vector with the computed ranks of the separable function.
        //! \details The operator \a _op should take three parameters: the first
        //!          is a reference to the ouput, the second the coefficient (or
        //!          sum of coefficient ) for that configuration element.
        //! \see CE::Mapping::VectorPlus::apply() for example.
        template< class T_COEFS, class T_VECIN, class T_VECOUT, class T_OP >
          void static rank_vector( const T_COEFS &_coefs, const T_VECIN &_vecin,
                                   T_VECOUT &_vecout, T_OP _op );
        //! Applies an operation throughout the separable function.
        template< class T_COEFS, class T_VECIN, class T_OUT, class T_OP >
          void static apply_throughout( const T_COEFS &_coefs, const T_VECIN &_vecin,
                                        T_OUT &_out, T_OP _op );
        //! \brief Applies an operation to dimension \a _d of the separable function.
        //! \details Output is a vector with each rank in a separate element.
        template< class T_COEFS, class T_VECIN, class T_VECOUT, class T_OP >
          void static apply_to_dim( const T_COEFS &_coefs, const T_VECIN &_vecin,
                                    T_VECOUT &_out, size_t _d, T_OP _op );
        //! \brief Applies an operation to rank \a _r of the separable function.
        //! \details Output is a scalar.
        template< class T_COEFS, class T_VECIN, class T_OUT, class T_OP >
          void static apply_to_rank( const T_COEFS &_coefs, const T_VECIN &_vecin,
                                     T_OUT &_out, size_t _r, T_OP _op );
        //! \brief Applies an operation to dimension \a _d and rank \a _r of the
        //!        separable function.
        //! \details Output is a scalar.
        template< class T_COEFS, class T_VECIN, class T_OUT, class T_OP >
          void static apply_to_dim_n_rank( const T_COEFS &_coefs, const T_VECIN &_vecin,
                                           T_OUT &_out, size_t _d, size_t _r, T_OP _op );
        //! Applies a normalization functor.
        template< class T_COEFS, class T_NORMS > 
          void static normalize( T_COEFS &_coefs, T_NORMS &_norms );
        //! Applies a randomization functor.
        template< class T_COEFS > 
          void static randomize( T_COEFS &_coefs, typename T_COEFS::value_type _howrandom );
      };

      //! Interface between the coefficient matrix and the separable function.
      class DirectCoefficients : public Indirection::Object
                                        < boost::numeric::ublas::matrix<types::t_real> >
      {
        typedef Indirection::Object
                < boost::numeric::ublas::matrix<types::t_real> > t_Base;
        public:
          //! Type of the matrix.
          typedef t_Base :: t_Object t_Matrix;
   
          //! Resizes matrix.
          void resize( size_t _i, size_t _j )
            { this->operator()().resize(_i, _j); }
      };
      
      //! Interface between the coefficient matrix and the separable function.
      class MatrixRangeCoefficients
      {
        //! Type of the matrix.
        typedef boost::numeric::ublas::matrix<types::t_real> t_RealMatrix;

        public:
          //! Type of the coefficients.
          typedef boost::numeric::ublas::matrix_range< t_RealMatrix > t_Matrix;
   
          //! Constructor.
          MatrixRangeCoefficients() {}
          //! Copy Constructor.
          MatrixRangeCoefficients   ( const MatrixRangeCoefficients &_c )
                                  : matrix_(_c.matrix_) {} 

          //! Returns a reference to the matrix.
          t_Matrix& operator()() { return *matrix_; }
          //! Returns a reference to the matrix.
          const t_Matrix& operator()() const { return *matrix_; }

          //! Sets range and reference matrix.
          void set( t_RealMatrix& _mat,
                    const boost::numeric::ublas::range &_a,
                    const boost::numeric::ublas::range &_b )
          {
            if( not matrix_.get() )
              matrix_.reset( new t_Matrix( _mat, _a, _b ) ); 
            else 
            {
              boost::shared_ptr<t_Matrix> y( new t_Matrix( _mat, _a, _b ) );
              matrix_ = y; 
            }
          }
          //! Resize matrix.
          void resize( size_t, size_t ) {};
   
        protected:
          //! The coefficients themselves.
          boost::shared_ptr<t_Matrix> matrix_;
      };
    
    } // end of Policy namespace.
  }
} // namespace LaDa
#include "functional.impl.h"

#endif
