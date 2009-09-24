//
//  Version: $Id$
//
#ifndef LADA_ATOMIC_POTENTIAL_VARIABLE_SET_H_
#define LADA_ATOMIC_POTENTIAL_VARIABLE_SET_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/io.hpp>
#include <string>

#include "../numeric_types.h"
#include "../functions.h"
#include "fitting_set.h"
#include "values.h"

namespace LaDa
{
  // Forward declaration
  //! \cond
  namespace Crystal
  {
    template<class T_TYPE> class TStructure;
  }
  //! \endcond

# ifdef __DOPYTHON
  namespace Python
  {
    void expose_collapse();
  }
# endif 
 
  namespace atomic_potential
  {
    // Forward declaration.
    class SumOfSeparables;

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
      //!             P_{r,i}[x^(u,v)}_{j\neq i}] = \prod_{j\neq i} \sum_o
      //!               \alpha_{r,j,o} g_{r,j,o}( x_j^{(u,v)} )
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
#       ifdef __DOPYTHON 
          friend void Python::expose_collapse();
#       endif
        protected:
          //! Type of the container of coefficients.
          typedef std::vector<vector_type>  t_Coefficients;
          //! Type of container of scaling factors.
          typedef std::vector<numeric_type> t_ScalingFactors;
          //! Values associated with the fitting structures.
          typedef FittingSet t_FittingSet;
          //! Values associated with the fitting structures and functions of the sum of seps.
          typedef Values t_Values;

        public:
          //! Constructor.
          Collapse(SumOfSeparables & _sumofseps);
          //! Copy Constructor.
          Collapse   (Collapse const &_c)
                   : fitting_set_(_c.fitting_set_),
                     values_(_c.values_),
                     coefficients_(_c.coefficients_),
                     scaling_factors_(_c.scaling_factors_),
                     sumofseps_(_c.sumofseps_) {}
     
          //! Creates A matrix and b vector for coordinate \a _i. 
          template<class T_MATRIX, class T_VECTOR>
            void lsq_data(T_MATRIX& _matrix, T_VECTOR& _vector, size_t _i) const;
          //! Updates coordinate \a _i.
          void update(size_t _i) { values_.update(coefficients_[_i], fitting_set_, _i); }
          //! Adds a structure to the fitting set.
          void add(Crystal::TStructure<std::string> const &_structure);
          //! Reassigns coefficients.
          void reassign() const;
          //! Returns coefficients for coordinate \a _i.
          t_Coefficients::value_type& coefficients(size_t _i) { return coefficients_[_i]; }
          //! Returns coefficients for coordinate \a _i.
          t_Coefficients::value_type const& coefficients(size_t _i) const 
            { return coefficients_[_i]; }
          //! Returns the number of coordinates.
          size_t nb_coordinates() const;
          //! Returns average error per structure in fitting set.
          numeric_type convergence() const;
     
        private:
          //! Weight of each structure. 
          t_FittingSet fitting_set_;
          //! Weight of each structure. 
          t_Values values_;
          //! Coefficients of the separable function.
          t_Coefficients coefficients_;
          //! Coefficients of the separable function.
          t_ScalingFactors scaling_factors_;
          //! Sum of separables function.
          SumOfSeparables &sumofseps_;
      };



      template< class T_MATRIX, class T_VECTOR>
        void Collapse::lsq_data(T_MATRIX& _matrix, T_VECTOR& _vector, size_t _i) const
        {
          LADA_ASSERT( _i < coefficients_.size(),
                       "index out-of-range: " << _i << " >= " <<  coefficients_.size() << ".\n")
          const size_t vec_size( coefficients_[_i].size() );
          if( _matrix.size2() != vec_size )
          {
            _matrix.resize(vec_size, vec_size);
            _vector.resize(vec_size);
          }
          else if( _matrix.size1() != vec_size )
            _matrix.resize(vec_size, vec_size);
          for( size_t i(0); i < vec_size; ++i )
          {
            for( size_t j(0); j < vec_size; ++j )
              _matrix(i,j) = 0e0;
            _vector(i) = 0e0;
          }
          
          t_FittingSet :: str_iterator i_str = fitting_set_.begin(_i);
          t_FittingSet :: str_iterator const i_str_end = fitting_set_.end(_i);
          t_Values::const_str_iterator i_str_val = values_.begin(_i);
          typename T_VECTOR :: value_type const nb_structures( 1e0 / fitting_set_.size() );
          //! Loop over structures.
          for(; i_str != i_str_end; ++i_str, ++i_str_val)
          {
            numeric_type const str_weight(i_str.weight());
            if( str_weight == 0e0 ) continue; // don't fit this structure.
        
            t_FittingSet::str_iterator::rep_iterator i_rep = i_str.begin();
            t_FittingSet::str_iterator::rep_iterator const i_rep_end = i_str.end();
            t_Values::const_str_iterator::const_rep_iterator i_rep_val = i_str_val.begin();
            vector_type G(vec_size, 0); // values to sum to matrix and vector.
        
            // loop over structure representations.
            for(; i_rep != i_rep_end; ++i_rep, ++i_rep_val )
            {
              numeric_type const rep_weight(i_rep.weight() * str_weight);
        
              typedef t_Values::const_str_iterator::const_rep_iterator
                              ::const_rank_iterator rank_iterator;
              rank_iterator i_rank_val = i_rep_val.begin();
              rank_iterator const i_rank_val_end = i_rep_val.end();
              t_ScalingFactors::const_iterator i_scale = scaling_factors_.begin();
#             ifdef LADA_DEBUG
                t_ScalingFactors::const_iterator const i_scale_end = scaling_factors_.end();
#             endif
        
              // loop over ranks.
              for(size_t i(0); i_rank_val != i_rank_val_end; ++i_rank_val, ++i_scale )
              {
                LADA_ASSERT(i_scale != i_scale_end, "Iterator out of range.\n")
                numeric_type factor_i( i_rank_val.other() * (*i_scale) * rep_weight );
                rank_iterator::function_iterator i_func = i_rank_val.begin();
                rank_iterator::function_iterator const i_func_end = i_rank_val.end();
                // loop over inner functions.
                for(; i_func != i_func_end; ++i_func, i+=Functions::N )
                  G( i + i_rep.specie() ) += (*i_func) * factor_i;
              } // loop over ranks
            } // loop over representations
        
            // now adds to A matrix and b vector.
            for(size_t i(0); i < vec_size; ++i )
            {
              for(size_t j(i); j < vec_size; ++j)
                _matrix(i,j) += G(i) * G(j) * str_weight * nb_structures;
              _vector(i) += i_str.energy() * G(i) * str_weight * nb_structures;
            } 
          } // loop over structures.
        
          // creates second half of matrix.
          for( size_t i(1); i < vec_size; ++i )
            for( size_t j(0); j < i; ++j )
              _matrix(i,j) = _matrix(j,i);
        }
    } // namespace collapse
  } // namespace atomic_potential
} // namespace LaDa
#endif
