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
  namespace atomic_potential
  {
    namespace collapse
    {

      bool Collapse::lsq(matrix_type &_matrix, vector_type &_vector, size_t _i) const
      {
        namespace bnu = boost::numeric::ublas;

        const size_t vec_size( coefficients_[_i].size() );
        if( _matrix.size2() != vec_size )
        {
          _matrix.resize(vec_size, vec_size);
          _vector.resize(vec_size);
        }
        else if( _matrix.size2() != _vector.size() )
          _matrix.resize(vec_size, vec_size);
        for( size_t i(0); i < vec_size; ++i )
        {
          for( size_t j(0); j < vec_size; ++j )
            _matrix(i,j) = 0e0;
          _vector(i) = 0e0;
        }
          
        
        t_FittingStructures :: str_iterator i_str = fitting_structures_.begin(_i);
        t_FittingStructures :: str_iterator const i_str_end = fitting_structures_.end(_i);
        t_Values::str_iterator i_str_val = values_.begin(_i);
        //! Loop over structures.
        for(; i_str != i_str_end; ++i_str, ++i_str_val)
        {
          numeric_type const str_weight(i_str->weight());
          if( str_weight == 0e0 ) continue; // don't fit this structure.

          t_Representations::str_iterator::rep_iterator i_rep = i_str->begin();
          t_Representations::str_iterator::rep_iterator const i_rep_end = i_str->end();
          t_Values::str_iterator::rep_iterator i_rep_val = i_str_val->begin();
          vector_type G(vec_size, 0); // values to sum to matrix and vector.

          // loop over structure representations.
          for(; i_rep != i_rep_end; ++i_rep, ++i_rep_val )
          {
            numeric_type const rep_weight(i_rep->weight() * str_weight);

            typedef t_FittingStructures::str_iterator::rep_iterator
                                       ::coordinate_iterator typedef coordinate_iterator;
            typedef t_Values::str_iterator::rep_iterator
                            ::rank_iterator typedef rank_iterator;
            rank_iterator i_rank_val = i_rep_val->begin();
            rank_iterator const i_rank_val_end = i_rep_val->end();
            coordinate_iterator i_coord = i_rep->begin();

            // loop over ranks.
            for(size_t i(0); i_rank_val != i_rank_val_end; ++i_rank_val )
            {
              numeric_type factor_i( i_rank_val->other() );
              rank_iterator::function_iterator i_func = i_rank_val->begin();
              rank_iterator::function_iterator const i_func_end = i_rank_val->end();
              // loop over inner functions.
              for(; i_func != i_func_end; ++i_func, i+=Functions::N, ++i_coord)
                G( i + (*i_coord) ) += rep_weight * (*i_rank_val) * other;
            } // loop over ranks
          } // loop over representations

          // now adds to A matrix and b vector.
          for(size_t i(0); i < vec_size; ++i )
          {
            for(size_t j(i); j < vec_size; ++j)
              _matrix(i,j) += G(i) * G(j) * str_weight;
            _vector(i) += i_str->energy() * G(i) * str_weight;
          } 
        } // loop over structures.

        // creates second half of matrix.
        for( size_t i(0); i < vec_size; ++i )
          for( size_t j(i+1); j < vec_size; ++j )
            _matrix(i,j) = _matrix(j,i);
      }

    } // namespace collapse
  } // namespace atomic_potential
} // namespace LaDa
