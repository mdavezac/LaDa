//
//  Version: $Id$
//
#ifndef _SEPARABLE_COLLAPSE_H_
#define _SEPARABLE_COLLAPSE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<vector>
#include<functional>
#include<numeric>

#include <boost/lambda/lambda.hpp>
#include <boost/type_traits/function_traits.hpp>

#include<opt/types.h>
#include<opt/debug.h>

//! Contains all things separable functions.
namespace Separable
{
  /** \brief Collapses a sum of separable function into Fitting::Allsq format.
   *  \details This flavor keeps track of computed basis functions
   *           (\f$g_{d,i}^{(r)}\f$).  This is quite a complex process. A sum
   *           of separable functions is expressed as follows: 
   *    \f[
   *        F(\mathbf{x}) = \sum_r \prod_d \sum_i \lambda_{d,i}^{(r)}
   *                        g_{i,n}^{(r)}(x_i),
   *    \f]
   *    We will keep track of the r, i, and n indices. Furthermore, the least
   *    square-fit method fits to  \e o targets values. It expects a "2d"
   *    input matrix I[d, (r,i) ], where the slowest running index is the
   *    dimension \e d.  I[d] should be a vecrtor container type of scalar
   *    values. These scalar values are orderer with  \e i the fastest running
   *    index and \e r the slowest. The method also expects a functor which
   *    can create the matrix \e A (as in \e Ax = \e b ) for each specific
   *    dimension \e d. The type of A is a vector of scalars. The fastest
   *    running index is \e i, followed by \e r, and finally \e o. 
   *    \tparam T_FUNCTION must be a Function< Summand< T_BASIS > > 
   **/
  template< class T_FUNCTION >
  class Collapse
  {
    protected:
      //! Return type of the function
      typedef typename T_FUNCTION :: t_Return t_Type;

    public:
      //! Type of the sum of separable functions to collapse.
      typedef T_FUNCTION t_Function;
      //! Reference to the function to fit.
      T_FUNCTION &function;

      //! Constructor
      Collapse   ( t_Function &_function )
                  : D(0), nb_targets(0), nb_ranks(0), function( _function ) {}

      /** \brief Constructs the matrix \a A for dimension \a _dim. 
       *  \tparam T_VECTOR is a vector.
       *  \tparam T_MATRIX is a vector, e.~g. a flattened matrix. This class is
       *                   meant for %GSL which expects this type of memory
       *                   setup for matrices. 
       *  \tparam T_VECTORS is a vector of vectors or similar. 
       *  \param[out] _b are the target value of the collapsed least square
       *                 fit along dimension \a _dim, \a _b[ (r,i) ].
       *      \f[
       *         \_b[ (r,i) ] = \sum_o \text{targets}[o] g_{i,d}^{(r)}(x_d^{(o)})
       *                              \prod_{u\neq d} \text{factor}[u,(r,o)]
       *      \f]
       *  \param[out] _A will contain the collapsed matrix for the linear
       *                 least-square fit. 
       *      \f[
       *         \_A[ (r,i,r',i') ] = \sum_o g_{i,d}^{(r)}(x_d^{(o)})
       *                              \prod_{u\neq d} \text{factors}[u,(r,o)]
       *                              \prod_{u\neq d} \text{factors}[u,(r',o)]
       *                               g_{i',d}^{(r')}(x_d^{(o)})
       *      \f]
       *  \param[in] The dimension for which to create \a _A.
       *  \param[in] _targets are the structural target value, \a _target[o].
       *  \param[in] _coefs contains the coefficients for \e all dimensions.
       *                    \a _coefs[d, (r,i) ]
       **/
      template< class T_VECTOR, class T_MATRIX, class T_VECTORS >
      void operator()( T_VECTOR &_b, T_MATRIX &_A, types::t_unsigned _dim,
                       const T_VECTOR &_targets, const T_VECTORS &_coefs  );

      //! \brief Constructs the completely expanded matrix.
      //! \tparams T_VECTORS is a vector of vectors.
      //! \tparams T_VECTOR is a vector of return types.
      //! \params[in] _x contains the input to sum of separable function in the
      //!                format \a _x[ o,d ].
      //! \params[in] _w contains the weights attached to each structural
      //!                target value. This array is copied to a local
      //!                container. \a _w[o].
      template< class T_VECTORS, class T_VECTOR > 
        void init( const T_VECTORS &_x, const T_VECTOR &_w );
      //! Resets collapse functor. Clears memory.
      void reset();
      //! Creates a collection of random coefficients.
      //! \details A vector of vectors of correct dimensions is created and
      //!          initialized.
      //! \tparams should be a vector of vectors.
      //! \params[out] _coefs creates \a _coefs[d, (r,i) ] uninitialized.
      template< class T_VECTORS > void create_coefs( T_VECTORS &_coefs,
                                                     t_Type _howrandom= 0.5e0 ) const;
      //! Assigns solution coefficients to Collapse::function.
      //! \tparam T_VECTORS is a vector of vectors or similar. 
      //! \param[in] _solution contains the coefficients for \e all dimensions.
      //!                      \a _solution[d, (r,i) ]
      template< class T_VECTORS > void reassign( const T_VECTORS &_solution ) const;
      //! Evaluates the sum of squares.
      template< class T_VECTORS, class T_VECTOR >
        typename t_Function::t_Return evaluate( const T_VECTORS &_coefs,
                                                const T_VECTOR &_targets );

      //! \brief Initializes Collapse::factors using values from \a _coefs.
      //! \details Also normalizes coefficients (r,d) and stores norm in rank
      //!          coefficient.
      //! \tparam T_VECTORS is a vector of vectors or similar. 
      //! \param[in] _coefs contains the coefficients for \e all dimensions.
      //!                   \a _coefs[d, (r,i) ]
      template< class T_VECTORS > void update_all( T_VECTORS &_coefs );
      //! \brief Updates Collapse::factors for dimension \a _dim.
      //! \details Also normalizes coefficients (r,d) and stores norm in rank
      //!          coefficient.
      //! \tparam T_VECTORS is a vector of vectors or similar. 
      //! \param[in] The dimension for which to create \a _A. Passed from
      //!            Collapse::operator()().
      //! \param[in] _coefs contains the coefficients for \e all dimensions.
      //!                   \a _coefs[d, (r,i) ].
      template< class T_VECTORS >
        void update( types::t_unsigned _dim, T_VECTORS &_coefs );

    protected:
      //! Maximum dimension.
      types::t_unsigned D;
      //! Number of target values.
      types::t_unsigned nb_targets;
      //! Number of ranks.
      types::t_unsigned nb_ranks;
      //! \brief Type of the matrix containing expanded function elements.
      typedef std::vector< std::vector< std::vector< t_Type > > > t_Expanded;
      //! A matrix with all expanded function elements.
      //! \details expanded[ d, r, (i,o) ]. The type is a vector of vectors.
      //!          The fastest-running \e internal index is \e i. The
      //!          fastest-running \e external index is d.
      //!          \f$ \text{expanded}[ d, r, (i,o)] = g_{d,i}^{(r)}(x_d^{(o)})\f$.
      t_Expanded expanded;
      //! Type of the factors.
      typedef std::vector< std::vector< t_Type > > t_Factors;
      /** \brief A matrix wich contains the factors of the separable functions.
       *  \details factors[ (r, o), d ]. A vector of vectors. 
       *           \f$
       *               \text{factors}[(r, o), d] = \sum_i
       *               \lambda_{d,i}^{(r)}g_{d,i}^{(r)}(x_d^{(o)})
       *           \f$. **/
      t_Factors factors;
      //! Type of the sizes.
      typedef std::vector< std::vector< types::t_unsigned > > t_Sizes;
      //! \brief Sizes of the basis per dimension and rank.
      //! \details sizes[ d, r ] = max i \@ (d,r). r is the fastest running
      //!          index.
      t_Sizes sizes;
      //! Type of the weights.
      typedef std::vector< t_Type > t_Weights;
      //! Weights attached to each structural target value.
      t_Weights weights;
  };

} // end of Separable namespace

#include "collapse.impl.h"

#endif
