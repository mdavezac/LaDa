//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/vector_proxy.hpp>

#include "collapse.h"

namespace LaDa
{
  namespace AtomicPotential
  {
    bool Collapse::lsq(t_Matrix &_matrix, t_Vector &_b, size_t _i) const
    {
      namespace bnu = boost::numeric::ublas;
      // Coefficients of this variable.
      t_SeparableSizes::const_iterator i_size( separable_size.begin() );
      t_SeparableSizes::const_iterator const i_size_end( separable_size.end() );
      types::t_real const function_size( variable_major_.coefficients_[_i].size() );
      
      //! Loop over structures.
      t_Weights :: const_iterator i_weight = weights.begin();
      t_Weights :: const_iterator const i_weight_end = weights.end();
      t_SymWeights :: const_iterator i_sym_weight;
      t_SymWeights :: const_iterator const i_sym_weight_end = weights.end();
      t_FunctionValues::value_type::const_iterator i_str_val = function_values_[_i].begin();
      t_Factors::value_type::const_iterator i_str_factor( factors_[_i].begin() ); 
      t_Values::value_type::const_iterator i_str_values( values_[_i].begin() ); 
      for(; i_weight != i_weight_end; ++i_weight, ++i_str_val, ++i_str_factor) // loop over structures.
      {
        t_Vector G(function_size, 0);
        t_FunctionValues::value_type::value_type::const_iterator i_sym_val = i_str_val->begin();
        t_Factors::value_type::value_type::const_iterator i_sym_factor( i_str_factor->begin() );
        t_Values::value_type::value_type::const_iterator i_sym_values( i_str_values->begin() );
        i_sym_weight = weights.begin();
        // loop over symmetrized structures.
        for(; i_sym_weight != i_sym_weight_end; ++i_sym_weight, ++i_sym_val, ++i_str_factor)
        {
          t_SeparableSizes::const_iterator i_size = separable_sizes_[_i].begin();
          t_SeparableSizes::const_iterator i_size_end = separable_sizes_[_i].end();
          for(size_t start(0); i_size != i_size_end; ++i_size) // loop over ranks.
          {
            size_t const end( start + (*i_size) );
            bnu::vector_range<t_Vector> x(left, bnu::range(start, end));
            bnu::vector_range<t_Vector> coefs(, bnu::range(start, end));
            range *= 
            
            start = *i_size;
          }
        } // loop over symmetrized structures.
      } // loop over structures.
    }

  } // namespace AtomicPotential
} // namespace LaDa
