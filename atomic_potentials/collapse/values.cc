//
//  Version: $Id: collapse.cc 1270 2009-08-17 03:11:58Z davezac $
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "../sum_of_separables.h"
#include "../representation.h"
#include "values.h"
#include "fitting_set.h"

namespace LaDa
{
  namespace atomic_potential
  {
    namespace collapse
    {
      Values::str_iterator Values::begin(size_t _i) 
      {
        LADA_ASSERT( _i < coord_rank_values_.size(), "Index out of range.\n" );
 
        str_iterator result(_i, coord_rank_values_);
        result.i_coord_rank_value_ = coord_rank_values_[_i].begin();
        result.i_rank_value_ = rank_values_.begin();
        result.i_function_values_ = function_values_[_i].begin();
        return result;
      }
      Values::str_iterator Values::end(size_t _i) 
      {
        LADA_ASSERT( _i < coord_rank_values_.size(), "Index out of range.\n" );
 
        str_iterator result(_i, coord_rank_values_);
        result.i_coord_rank_value_ = coord_rank_values_[_i].end();
        // result.i_rank_value_ = rank_values_.end();
        // result.i_function_values_ = function_values_[_i].end();
        return result;
      }
      Values::const_str_iterator Values::begin(size_t _i) const
      {
        LADA_ASSERT( _i < coord_rank_values_.size(), "Index out of range.\n" );
 
        const_str_iterator result(_i, coord_rank_values_);
        result.i_coord_rank_value_ = coord_rank_values_[_i].begin();
        result.i_rank_value_ = rank_values_.begin();
        result.i_function_values_ = function_values_[_i].begin();
        return result;
      }
      Values::const_str_iterator Values::end(size_t _i) const
      {
        LADA_ASSERT( _i < coord_rank_values_.size(), "Index out of range.\n" );
 
        const_str_iterator result(_i, coord_rank_values_);
        result.i_coord_rank_value_ = coord_rank_values_[_i].end();
        // result.i_rank_value_ = rank_values_.end();
        // result.i_function_values_ = function_values_[_i].end();
        return result;
      }

      void Values::update(vector_type const& _coefs, t_FittingSet const &_ftstr, size_t _i )
      {
        typedef SumOfSeparables::const_coord_range const_coord_range;
        typedef const_coord_range::const_rank_range const_rank_range;
        typedef const_rank_range::const_iterator const_iterator_functions;
        t_FittingSet :: str_iterator i_str = _ftstr.begin(_i);
        t_FittingSet :: str_iterator const i_str_end = _ftstr.end(_i);
        //! Loop over structures.
        for(size_t s(0); i_str != i_str_end; ++i_str, ++s)
        {
          if( i_str.weight() == 0e0 ) continue; // don't fit this structure.
      
          t_FittingSet::str_iterator::rep_iterator i_rep = i_str.begin();
          t_FittingSet::str_iterator::rep_iterator const i_rep_end = i_str.end();
          // loop over structure representations.
          for(size_t p(0); i_rep != i_rep_end; ++i_rep, ++p )
          {
      
            // loop over ranks.
            size_t const Nrank(rank_values_[s][p].size());
            for(size_t r(0), c(0); r < Nrank; ++r )
            {
              numeric_type &coord_rank_value(coord_rank_values_[_i][s][p][r]);
              numeric_type &rank_value(rank_values_[s][p][r]);
              types::t_real i_neq_j(1);
              if( Fuzzy::is_zero(coord_rank_value) )
                for(size_t i(0), Nvar(coord_rank_values_.size()); i < Nvar; ++i)
                  if( i != _i ) i_neq_j *= coord_rank_values_[i][s][p][r];
              else i_neq_j = rank_value / coord_rank_value;
      
              // Updates all i dependent factors.
              coord_rank_value = numeric_type(0);
              typedef t_FunctionValues::value_type::value_type
                                      ::value_type::value_type::const_iterator t_Functions;
              t_Functions i_func( function_values_[_i][s][p][r].begin() );
              t_Functions const i_func_end( function_values_[_i][s][p][r].end() );
              for(; i_func != i_func_end; ++i_func, c+=Functions::N )
                coord_rank_value += _coefs(c + i_rep.specie()) * (*i_func);
      
              rank_value = i_neq_j * coord_rank_value;
            } // loop over ranks
          } // loop over representations
        } // loop over structures.
      }


      void Values::add( Representation const &_reps,
                        SumOfSeparables const& _sumofseps ) 
      {
        typedef SumOfSeparables::const_coord_range const_coord_range;
        typedef const_coord_range::const_rank_range const_rank_range;
        typedef const_rank_range::const_iterator const_iterator_functions;
        const_coord_range const _range(_sumofseps.range());

        // Case when first structure to be added to function.
        if( function_values_.size() == 0 ) function_values_.resize( _range.size() );
        if( coord_rank_values_.size() == 0 ) coord_rank_values_.resize( _range.size() );
          
        // Sanity checks.
        LADA_ASSERT( _range.size() == function_values_.size(),
                     "Inconsistent container size.\n" )
        LADA_ASSERT( _range.size() == coord_rank_values_.size(),
                     "Inconsistent container size.\n" )

        { // resizes containers.
          t_FunctionValues::iterator i_funcs = function_values_.begin();
          t_FunctionValues::iterator const i_funcs_end = function_values_.end();
          t_CoordRankValues::iterator i_coords = coord_rank_values_.begin();
          for(; i_funcs != i_funcs_end; ++i_coords, ++i_funcs)
          {
            i_funcs->resize( i_funcs->size()+1 );
            i_coords->resize( i_coords->size()+1 );
          } // loop over variables.
          rank_values_.resize(rank_values_.size() + 1);
        }
        t_RankValues::value_type &str_rv = rank_values_.back();

        // Loops over representations.
        Representation::const_iterator i_rep = _reps.begin();
        Representation::const_iterator const i_rep_end = _reps.end();
        for(; i_rep != i_rep_end; ++i_rep)
        {
          // loop overvariables
          VariableSet::t_Variables::const_iterator i_var = i_rep->variables.begin(); 
#         ifdef LADA_DEBUG
            VariableSet::t_Variables::const_iterator const
              i_var_end = i_rep->variables.end();
#         endif
          for(const_coord_range range(_range); range; ++range, ++i_var)
          {
            LADA_ASSERT( i_var != i_var_end, "Iterator out of range.\n" )
            LADA_ASSERT( *range < coord_rank_values_.size(), "Index out-of-range.\n")
            LADA_ASSERT( *range < function_values_.size(), "Index out-of-range.\n")
            t_CoordRankValues::value_type::value_type
               &rep_crv( coord_rank_values_[*range].back() );
            t_FunctionValues::value_type::value_type
               &rep_fv( function_values_[*range].back() );
            
            // Adds rep to coord and str. 
            rep_crv.resize(rep_crv.size() + 1);

            t_FunctionValues::value_type::value_type::value_type funcs_rank;
            for(const_rank_range rank(range.range()); rank; ++rank)
            {
              t_FunctionValues::value_type::value_type
                              ::value_type::value_type funcs;
              const_iterator_functions i_func = rank.begin(); 
              const_iterator_functions i_func_end = rank.end(); 
              numeric_type var_value(0);
              for(; i_func != i_func_end; ++i_func)
              {
                funcs.push_back( (*i_func)( i_var->first ) );
                var_value += (*i_func)(*i_var);
              }
              funcs_rank.push_back(funcs);
              rep_crv.back().push_back(var_value);
            } // over ranks.
            rep_fv.push_back(funcs_rank);
          } // over variables.
          
          str_rv.resize(str_rv.size() + 1);
          t_RankValues::value_type::value_type &rep_rv = str_rv.back();
          SumOfSeparables::const_iterator i_sep = _sumofseps.begin();
          SumOfSeparables::const_iterator const i_sep_end = _sumofseps.end();
          for(; i_sep != i_sep_end; ++i_sep)
            rep_rv.push_back( i_sep->get<0>()(i_rep->variables) );
        } // over representations.
      } // end of add


    } // namespace collapse
  } // namespace atomic_potential
} // namespace LaDa
