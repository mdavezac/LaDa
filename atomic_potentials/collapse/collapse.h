//
//  Version: $Id$
//
#ifndef LADA_ATOMIC_POTENTIAL_VARIABLE_SET_H_
#define LADA_ATOMIC_POTENTIAL_VARIABLE_SET_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <boost/numeric/ublas/matrix.hpp>

#include "../numeric_types.h"

namespace LaDa
{
  namespace atomic_potential
  {
    // Forward declaration.
    class SumOfSeparables;
    class FittingSet;
    class Values;

    namespace collapse
    {
      //! \brief Linear least-square-variables input variables.
      //! \details The separable function is defined as:
      //!          \f[
      //!              S(X) = \sum_r s_r \prod_i \sum_o \alpha_{r,i,o} g_{r,i,o}( x_i )
      //!          \f]
      //!          With \a u indexing the input structures, and \a v indexing the
      //!          sum of symmetric structures, a least-square fit for variable
      //!          \a _i can be performed using the vector:
      //!          \f[
      //!             G^{(u)}_i=\sum_v w_v^{(u)} P_{r,i}[x_i^{(u,v)}] g_{i,r,o}[x_i^{(u,v)}]
      //!          \f],
      //!          where 
      //!          \f[
      //!             P_{r,i}[x^(u,v)}_{j\neq i}] = \prod_{j\neq i} \sum_o \alpha_{r,j,o} g_{r,j,o}( x_j^{(u,v)} )
      //!          \f].
      //!          To help with this process, containers over the fitting
      //!          structures and over function values are created,
      //!          with a hierarchy of iterators. These iterators go from loops
      //!          over fitting structures, to loops over representations of a
      //!          single structure, to loop over ranks, to a loop over inner
      //!          functions. They give access to the values required for
      //!          computing \f$G${(u)}\f$. 
      class Collapse
      {
        protected:
          //! Type of the container of coefficients.
          typedef std::vector<vector_type>  t_Coefficients;
          //! Type of container of scaling factors.
          typedef std::vector<SumOfSeparables::result_type> t_ScalingFactors;
          //! Values associated with the fitting structures.
          typedef FittingSet t_FittingSet;
          //! Values associated with the fitting structures and functions of the sum of seps.
          typedef Values t_Values;
          //! Type of the variable major version of the separable function.
          typedef VariableMajor t_VariableMajor;

        public:
          //! Constructor.
          Collapse(SumOfSeparables const &_sumofseps) : variable_major_(_sumofseps) {}
          //! Copy Constructor.
          Collapse   (Collapse const &_sos)
                   : fitting_set_(_c.fitting_set_),
                     values_(_c.values_),
                     coefficients_(_c.coefficients_),
                     scaling_factors_(_c.scaling_factors_) {}
     
          //! Creates A matrix and b vector for coordinate \a _i. 
          bool lsq_data(matrix_type &_matrix, vector_type &_vector, size_t _i) const;
          //! Return coefficients for coordinate \a _i.
          t_Coefficients::value_type& coefficients(size_t _i) { return coefficients_[_i]; }
          //! Updates coordinate \a _i.
          void update(size_t _i) { values_.update(_coefficients[_i], fitting_set_, _i); }
          //! Adds a structure to the fitting set.
          void add(Crystal::TStructure<std::string> const &_structure);
     
        private:
          //! Weight of each structure. 
          t_FittingSet fitting_set_;
          //! Weight of each structure. 
          t_Values values_;
          //! Coefficients of the separable function.
          t_Coefficients coefficients_;
          //! Coefficients of the separable function.
          t_ScalingFactors scaling_factors_;
          //! Variable major version of the separable function.
          t_VariableMajor variable_major_;
      };
    } // namespace collapse
  } // namespace atomic_potential
} // namespace LaDa
#endif
