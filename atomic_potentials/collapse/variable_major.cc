//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <iterator>

#include "../sum_of_separables.h"
#include "variable_major.h"

namespace LaDa
{
  namespace AtomicPotential
  {
    VariableMajor::VariableMajor(SumOfSeparables const& _sumofseps) 
    {
      try
      {
        *this = _sumoseps;
      }
      catch(...)
      {
        std::cerr << "Could not create VariableMajor.\n"; 
        coefficients_.clear();
        functions_.clear();
        scales_.clear();
      }
    }

    bool VariableMajor::operator=(SumOfSeparables const& _sumofseps) 
    {
      coefficients_.clear();
      functions_.clear();
      scales_.clear();

      scales_ = _sumofseps.coefficients_;

      //! Gets maximum number of variables 
      size_t max_vars(0);
      SumOfSeparables::t_Functions::const_iterator i_seps( _sumofseps.functions_.begin() );
      SumOfSeparables::t_Functions::const_iterator const i_seps_end( _sumofseps.functions_.end() );
      for(; i_seps != i_seps_end; ++i_seps)
        max_vars = std::max(max_vars, i_seps->size());

      //! loops over variables.
      coefficients_.resize( max_vars);
      functions_.resize( max_vars);
      t_Coefficients :: iterator i_coef( coefficients_.begin() );
      t_Coefficients :: iterator const i_coef_end( coefficients_.end() );
      t_Functions :: iterator i_func( functions_.end() );
      for(size_t i(0); i_coef != i_coef_end; ++i_coef, ++i_func, ++i)
      {
        // loops over ranks.
        i_seps = _sumofseps.functions_.begin();
        for(; i_seps != i_seps_end; ++i_seps)
        {
          if( i_seps->size() < i ) continue;

          // finds ith variable.
          typedef SumOfSeparables::t_Functions::value_type t_Separable;
          t_Separable::const_iterator i_var(i_seps->functions_.begin());
          for( size_t j(0); j < i; ++j, ++i_var);

          // Adds coefficients and functions.
          std::copy( i_var->functions_.begin(), i_var->functions_.end(), 
                     std::back_inserter(*i_func) );
          for(size_t j(0); j < Functions::N; ++j)
            for(size_t o(0); o < i_var->coefficients.size(); ++o)
              i_coef->push_back( i_var->coefficients[j][o] );
        } // loop over ranks.
      } // loop over variables.
      return true;
    }

    bool VariableMajor::reassign(SumOfSeparables & _sumofseps) const
    {
      LADA_ASSERT( _sumofseps.coefficients_.size() == scales_.size(), "Incoherent containers.\n" )
      scales_ = _sumofseps.coefficients_;

      //! loops over variables.
      t_Coefficients :: const_iterator i_coef( coefficients_.begin() );
      t_Coefficients :: const_iterator const i_coef_end( coefficients_.end() );
      t_Functions :: const_iterator i_func( functions_.end() );
      for(size_t i(0); i_coef != i_coef_end; ++i_coef, ++i_func, ++i)
      {
        // loops over ranks.
        i_seps = _sumofseps.functions_.begin();
        size_t current_index(0);
        for(; i_seps != i_seps_end; ++i_seps)
        {
          if( i_seps->size() < i ) continue;

          // finds ith variable.
          typedef SumOfSeparables::t_Functions::value_type t_Separable;
          t_Separable::const_iterator i_var(i_seps->functions_.begin());
          for( size_t j(0); j < i; ++j, ++i_var);

          // Adds coefficients and functions.
          size_t const N( i_var->functions_.size() );
          LADA_ASSERT( current_index + N <= i_coefs->size(), "Incoherent containers.\n" );
          std::copy
          ( 
            i_func->begin() + current_index,
            i_func->end() + current_index + N,
            i_var->functions_.begin()
          );
          for(size_t u(0), o(0); u < i_coef->size(); ++o)
            for(size_t j(0); j < Functions::N; ++j, ++u)
            {
              LADA_ASSERT( u < i_coef->size, "Index out of range.\n" )
              i_var->coefficients[j][o] = (*i_coef)[u];
            }
          current_index += N;
        } // loop over ranks.
#       ifdef LADA_DEBUG
          for(size_t j(0); j < Functions::N; ++j)
            LADA_ASSERT( current_index == (*i_coef)[i].size(), "Incoherent containers.\n" )
#       endif
      } // loop over variables.
      return true;
    }

  } // namespace AtomicPotential
} // namespace LaDa
