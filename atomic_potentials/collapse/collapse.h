//
//  Version: $Id$
//
#ifndef LADA_ATOMIC_POTENTIAL_VARIABLE_SET_H_
#define LADA_ATOMIC_POTENTIAL_VARIABLE_SET_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/matrix.hpp>

#include "sum_of_separables.h"

namespace LaDa
{
  namespace atomic_potential
  {
    //! \brief Linear least-square-variables input variables.
    //! \details The separable function is defined as:
    //           \f[
    //               S(X) = \sum_r s_r \prod_i \sum_o \alpha_{r,i,o} g_{r,i,o}( x_i )
    //           \f]
    //           With \a u indexing the input structures, and \a v indexing the
    //           sum of symmetric structures, a least-square fit for variable
    //           \a _i can be performed using the vector:
    //           \f[
    //              G^{(u)}_i=\sum_v w_v^{(u)} P_{r,i}[x_i{(u,v)}] g_{i,r,o}[x_i^{(u,v)}]
    //           \f],
    //           where 
    //           \f[
    //              P_{r,i}[x_v^{(u)}] = \prod_{j\neq i} \sum_o \alpha_{r,j,o} g_{r,j,o}( x_j^{(u,v)} )
    //           \f],
    class Collapse
    {
        //! \brief Values of each function of a scalar in the sum of separable functions.
        //! \details Indices are arranged from outermost to innermost as:
        //!              _ i (over variables)
        //!              _ u (over structures)
        //!              _ v (over symmetrics of a structure)
        //!              _ r, o (over ranks and function of a scalar of each sum.
        //!              .
        typedef std::vector // over variables.
                <
                  std::vector // over input structures
                  <
                    std::vector // over symmetries of input structure.
                    <
                      SumOfSeparables::t_Coefficients // over rank and innermost sum.
                    >
                  >
                > t_FunctionValues;
        //! Sizes of the separable functions.
        typedef std::vector< std::vector<size_t> > t_SeparableSizes; 
        //! Weight of each structure.
        typedef std::vector<SumOfSeparables::result_type> t_Weights;
        //! Weight of each symmetrized structure.
        typedef std::vector< t_Weights > t_SymWeights;
        //! Value of each separable function.
        typedef std::vector<SumOfSeparables::result_type> t_SeparableValues;
      public:
        //! Type of the matrix argument.
        typedef boost::numeric::ublas::matrix<SumOfSeparables::result_type> t_Matrix;
        //! Type of the vector argument.
        typedef boost::numeric::ublas::matrix<SumOfSeparables::result_type> t_Vector;
        //! Constructor.
        Collapse(SumOfSeparables const &_sumofseps);
        //! Copy Constructor.
        Collapse   (Collapse const &_sos)
                 : function_values_(_c.function_values_),
                   separable_sizes_(_c.separable_sizes_),
                   weights_(_c.weights_),
                   sym_weights_(_c.sym_weights_),
                   separable_values(_c.separable_values) {}

        //! Creates A matrix for variable \a _i
        bool A(t_Matrix &_matrix, size_t _i) const;
        //! Creates b vector for variable \a _i
        bool b(t_Vector &_vector, size_t _i) const;

      private:
        //! \brief Values of each function for each possible input.
        //! \details \f$g_{r,i,o}(x_v^{(u)})\f$. 
        t_FunctionValues g_;
        //! Size of each separable function;
        t_SeparableSizes separable_sizes_;
        //! Weight of each structure. 
        t_Weights weights_;
        //! Weight of each inner structure. 
        t_SymWeights sym_weights_;
        //! Value of each separable function.
        t_SeparableValues separable_values;
        //! Variable major function.
        VariableMajor variable_major;
    };
  } // namespace atomic_potential
} // namespace LaDa
#endif
