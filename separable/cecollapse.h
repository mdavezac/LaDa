//
//  Version: $Id$
//
#ifndef _SEPARABLE_EQUIV_COLLAPSE_H_
#define _SEPARABLE_EQUIV_COLLAPSE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<vector>
#include<functional>
#include<numeric>

#include<boost/lambda/lambda.hpp>
#include<boost/type_traits/function_traits.hpp>

#include<opt/types.h>
#include<opt/debug.h>

#include "collapse.h"

//! Contains all things separable functions.
namespace Separable
{
  /** \brief Collapses a sum of separable function into Fitting::Allsq format,
   *         while keeping track of equivalent configurations.
   **/
  template< class T_FUNCTION >
  class EquivCollapse : public Collapse<T_FUNCTION>
  {
    protected:
      //! Type of the base class.
      typedef Collapse<T_FUNCTION> t_Base;
      //! Return type of the function
      typedef typename t_Base :: t_Type t_Type;
    public:
      //! Type of the separable function.
      typedef T_FUNCTION t_Function;
      using t_Base :: function;
#     ifdef __DOHALFHALF__
        //! A regularization factor.
        t_Type regular_factor;
        //! Constructor
        EquivCollapse   ( t_Function &_function )
                      : t_Base( _function ), regular_factor(0e0),
                        include_rank_coef_in_reg( true ) {}
#     else
        //! Constructor
        EquivCollapse ( t_Function &_function ) : t_Base( _function ) {}
#     endif

      //! Destructor
      ~EquivCollapse() {}

      //! \brief Constructs the matrix \a A for dimension \a _dim. 
      //! \see Collapse::operator()().
      template< class T_VECTOR, class T_MATRIX, class T_VECTORS >
      void operator()( T_VECTOR &_b, T_MATRIX &_A, types::t_unsigned _dim,
                       const T_VECTOR &_targets, const T_VECTORS &_coefs  );

      //! Constructs the completely expanded matrix.
      template< class T_VECTORS, class T_VECTOR, class T_EWEIGHTS > 
        void init( const T_VECTORS &_x, const T_VECTOR &_w,
                   const T_EWEIGHTS &_eweights );

      using t_Base :: reset;
      using t_Base :: update;
      using t_Base :: update_all;
#     ifdef __DOHALFHALF__
        //! Creates a collection of random coefficients.
        //! \details A vector of vectors of correct dimensions is created and
        //!          initialized.  This member is equivalent to the base class
        //!          member except when __DOHALFHALF__ is defined, in which
        //!          case, the True/True basis function is centered around 0
        //!          (rather than 1).
        //! \tparams should be a vector of vectors.
        //! \params[out] _coefs creates \a _coefs[d, (r,i) ] uninitialized.
        template< class T_VECTORS > void create_coefs( T_VECTORS &_coefs,
                                                       t_Type _howrandom= 0.5e0 ) const;
#     endif
      using t_Base :: reassign;
      //! Evaluates the sum of squares.
      template< class T_VECTORS, class T_VECTOR >
        typename t_Function::t_Return evaluate( const T_VECTORS &_coefs,
                                                const T_VECTOR &_targets );

    protected:

      using t_Base :: D;
      using t_Base :: nb_targets;
      using t_Base :: nb_ranks;
      //! \brief Type of the matrix containing expanded function elements.
      typedef typename t_Base::t_Expanded t_Expanded;
      using t_Base :: expanded;
      //! Type of the factors.
      typedef typename t_Base :: t_Factors t_Factors;
      using t_Base :: factors;
      //! Type of the sizes.
      typedef typename t_Base :: t_Sizes t_Sizes;
      using t_Base :: sizes;
      //! Type of the weights.
      typedef typename t_Base :: t_Weights t_Weights;
      using t_Base :: weights;
      //! Type of the equivalent structure indexing.
      typedef std::vector< std::vector<types::t_real> > t_eWeights;
      //! Equivalent structure indexing.
      t_eWeights eweights;
  };

} // end of Separable namespace

#include "cecollapse.impl.h"

#endif
