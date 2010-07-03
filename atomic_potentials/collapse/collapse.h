#ifndef LADA_ATOMIC_POTENTIAL_VARIABLE_SET_H_
#define LADA_ATOMIC_POTENTIAL_VARIABLE_SET_H_

#include "LaDaConfig.h"

#include <string>
#include <vector>

#include <boost/tuple/tuple.hpp>

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
      //!             P_{r,i}[x^{(u,v)}_{j\neq i}] = \prod_{j\neq i} \sum_o
      //!               \alpha_{r,j,o} g_{r,j,o}( x_j^{(u,v)} )
      //!          \f].
      //!          To help with this process, containers over the fitting
      //!          structures and over function values are created,
      //!          with a hierarchy of iterators. These iterators go from loops
      //!          over fitting structures, to loops over representations of a
      //!          single structure, to loop over ranks, to a loop over inner
      //!          functions. They give access to the values required for
      //!          computing \f$G{(u)}\f$. 
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
                     sumofseps_(_c.sumofseps_),
                     nb_funcs_(_c.nb_funcs_) {}
     
          //! \brief Creates A matrix and b vector for coordinate \a _i.
          //! \details \a _A and \a _b are the input matrices to a conjugate
          //!          gradient method.
          template<class T_MATRIX, class T_VECTOR>
            void lsq_data(T_MATRIX& _matrix, T_VECTOR& _vector, size_t _i) const;
          //! \brief Computes linear least square matrix. 
          //! \details \a _A and \a _b are the input matrices to a linear lest-square fit method.
          template<class T_MATRIX, class T_VECTOR>
            void lstsq_data(T_MATRIX &_matrix, T_VECTOR &_vector, size_t _i) const;

          //! Updates coordinate \a _i.
          void update(size_t _i, vector_type const &_x);
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
          boost::tuples::tuple<numeric_type, numeric_type, numeric_type> errors() const;
          //! Returns weight of structure \a _i.
          numeric_type get_weight( size_t _i ) { return fitting_set_.get_weight(_i); }
          //! Sets weight of structure \a _i. Returns old weight.
          numeric_type set_weight( size_t _i, numeric_type &_w )
            { return fitting_set_.set_weight(_i, _w); }
          //! Number of structures.
          size_t nb_structures() const { return fitting_set_.size(); }
          //! Weighted sum of the squared energies.
          numeric_type y_squared() const;
          //! Weighted sum of the squared energies.
          numeric_type sum_w() const;
     
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
          //! Number of functions per coordinate and rank.
          std::vector< std::vector<size_t> > nb_funcs_;
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
          numeric_type nstr(0);
          //! Loop over structures.
          for(; i_str != i_str_end; ++i_str, ++i_str_val)
          {
            numeric_type const str_weight(i_str.weight());
            if( str_weight == 0e0 ) continue; // don't fit this structure.
            nstr += str_weight;
        
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
                _matrix(i,j) += G(i) * G(j) * str_weight;
              _vector(i) += i_str.energy() * G(i) * str_weight;
            } 
          } // loop over structures.
        
          // creates second half of matrix.
          // divides by total structure weights.
          numeric_type const inv_nstr( numeric_type(1) / nstr );
          for( size_t i(0); i < vec_size; ++i )
          {
            for( size_t j(i+1); j < vec_size; ++j )
            {
              _matrix(i,j) *= inv_nstr;
              _matrix(j,i) = _matrix(i,j);
            }
            _vector(i) *= inv_nstr;
            _matrix(i,i) *= inv_nstr;
          }
        }



      template< class T_MATRIX, class T_VECTOR>
        void Collapse::lstsq_data(T_MATRIX& _matrix, T_VECTOR& _vector, size_t _i) const
        {
          LADA_ASSERT( _i < coefficients_.size(),
                       "index out-of-range: " << _i << " >= " <<  coefficients_.size() << ".\n")
          size_t n(0);
          t_FittingSet :: str_iterator i_str = fitting_set_.begin(_i);
          t_FittingSet :: str_iterator const i_str_end = fitting_set_.end(_i);
          for(; i_str != i_str_end; ++i_str)
            if( i_str.weight() != 0e0 ) ++n;

          const size_t vec_size( coefficients_[_i].size() );
          if( _matrix.size2() != n )
          {
            _matrix.resize(n, vec_size);
            _vector.resize(n);
          }
          else if( _matrix.size1() != vec_size )
            _matrix.resize(n, vec_size);
          for( size_t i(0); i < n; ++i )
          {
            for( size_t j(0); j < vec_size; ++j )
              _matrix(i,j) = 0e0;
            _vector(i) = 0e0;
          }
          
          t_Values::const_str_iterator i_str_val = values_.begin(_i);
          //! Loop over structures.
          i_str = fitting_set_.begin(_i);
          for(n = 0; i_str != i_str_end; ++i_str, ++i_str_val, ++n)
          {
            if( i_str.weight() == 0e0 ) continue; // don't fit this structure.
            numeric_type const str_weight( std::sqrt(i_str.weight()) );
        
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
                  _matrix(n, i + i_rep.specie() ) += (*i_func) * factor_i;
              } // loop over ranks
            } // loop over representations
        
            _vector(n) += i_str.energy() * str_weight;
          } // loop over structures.
        }

    } // namespace collapse
  } // namespace atomic_potential
} // namespace LaDa
#endif
