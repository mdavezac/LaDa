//
//  Version: $Id: collapse.cc 1270 2009-08-17 03:11:58Z davezac $
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "values.h"

namespace LaDa
{
  namespace atomic_potential
  {
    namespace collapse
    {
      Values::str_iterator Values::begin(size_t _i) const
      {
        LADA_ASSERT( _i < coordinates_.size(), "Index out of range.\n" );

        str_iterator result(_i, coord_rank_values_);
        result.i_coord_rank_value_ = coord_rank_values_[_i].begin();
        result.i_rank_value_ = rank_values_.begin();
        result.i_function_values_ = function_values_[_i].begin();
        return result;
      }
      Values::str_iterator Values::end(size_t _i) const
      {
        LADA_ASSERT( _i < coordinates_.size(), "Index out of range.\n" );

        str_iterator result(_i, coord_rank_values_);
        result.i_coord_rank_value_ = coord_rank_values_[_i].end();
        // result.i_rank_value_ = rank_values_.end();
        // result.i_function_values_ = function_values_[_i].end();
        return result;
      }

      template< class T_IT2, class T_IT3, class T_CRV >
        void update_impl( t_RankValues::const_iterator::const_iterator _it1,
                          t_RankValues::const_iterator::const_iterator _it1_end,
                          T_IT2 _it2, T_IT3 _it3, vector_type const &_coefs,
                          T_CRV const &_crv )
        {
          numeric_type const coord_rank_value( *i_coord_rank_value_ );
          numeric_type rank_value; 
          if( Fuzzy::is_zero(coord_rank_value) ) 
          {
            rank_value = 1e0;
            t_CoordRankValues::const_iterator i_cr = coord_rank_values_.begin();
            t_CoordRankValues::const_iterator const i_cr_end = coord_rank_values_.end();
            for( size_t i(0); i_cr != i_cr_end; ++i_cr, ++i)
              if( i != _i ) result *= (*i_cr)[n_str_][n_rep_][n_];
          }
        }
                          

      template<class T_IT1, class T_IT2, class T_IT3, class T_CRV> 
        void update_impl( T_IT1 _i1, T_IT1 const &_it1_end, 
                          T_IT2 _i2, T_IT3 _it3, vector_type const &_coefs,
                          T_CRV const &_crv )
        {
          for(; _it1 != _it1_end; ++_it1, ++_it2, ++_it3 )
            update_impl( _i1->begin(), _i1->end(), _i2->begin(), _i3->begin(), _coefs, _crv );
        }

      void Values::update( vector_type const& _coefs, size_t _i )
      {
        update_impl( rank_values_->begin(), rank_values->end(),
                     coord_rank_values_[_i]->begin(), function_values_[_i]->begin(),
                     _coefs, coord_rank_values_[_i] );
      }

      Values::str_iterator::rep_iterator Values::str_iterator::begin() const
      {
        rep_iterator result(i_, n_, coord_rank_values_);
        result.i_coord_rank_value_ = i_coord_rank_value_->begin();
        result.i_rank_value_ = i_rank_value_->begin();
        result.i_function_values_ = function_values_->begin();
        return result;
      }
      Values::str_iterator::rep_iterator Values::str_iterator::end() const
      {
        rep_iterator result(i_, n_, coord_rank_values_);
        result.i_coord_rank_value_ = i_coord_rank_value_->end();
        // result.i_rank_value_ = i_rank_value_->end();
        // result.i_function_values_ = function_values_->end();
        return result;
      }

      Values::str_iterator::rep_iterator::rank_iterator
        Values::str_iterator::rep_iterator::begin() const
        {
          rank_iterator result(i_, n_str_, n_, coord_rank_values_);
          result.i_coord_rank_value_ = i_coord_rank_value_->begin();
          result.i_rank_value_ = i_rank_value_->begin();
          result.i_function_values_ = function_values_->begin();
          result.i_coefficients_ = coefficients_[i_].begin();
          return result;
        }
      Values::str_iterator::rep_iterator::rank_iterator
        Values::str_iterator::rep_iterator::end() const
        {
          rank_iterator result(i_, n_str_, n_, coord_rank_values_);
          result.i_coord_rank_value_ = i_coord_rank_value_->end();
          // result.i_rank_value_ = i_rank_value_->end();
          // result.i_function_values_ = function_values_->end();
          // result.i_coefficients_ = coefficients_[i_].end();
          return result;
        }

      numeric_type Values::str_iterator::rep_iterator::rank_iterator::other() const
      {
        numeric_type const coord_rank_value( *i_coord_rank_value_ );
        if( Fuzzy::is_zero(coord_rank_value) ) 
        {
          numeric_type result(1);
          t_CoordRankValues::const_iterator i_cr = coord_rank_values_.begin();
          t_CoordRankValues::const_iterator const i_cr_end = coord_rank_values_.end();
          for( size_t i(0); i_cr != i_cr_end; ++i_cr, ++i)
            if( i != _i ) result *= (*i_cr)[n_str_][n_rep_][n_];
          return result;
        }
        return (*i_rank_value) / coord_rank_value;
      }
    } // namespace collapse
  } // namespace atomic_potential
} // namespace LaDa
