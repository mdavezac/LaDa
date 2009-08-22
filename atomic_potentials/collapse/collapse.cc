//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/vector_proxy.hpp>

#include <crystal/structure.h>

#include "../sum_of_separables.h"
#include "../functions.h"
#include "../representation.h" 
#include "collapse.h"
#include "values.h"
#include "fitting_set.h"

namespace LaDa
{
  namespace atomic_potential
  {
    namespace collapse
    {

      bool Collapse::lsq_data(matrix_type &_matrix, vector_type &_vector, size_t _i) const
      {
        namespace bnu = boost::numeric::ublas;

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
        t_Values::str_iterator i_str_val = values_.begin(_i);
        //! Loop over structures.
        for(; i_str != i_str_end; ++i_str, ++i_str_val)
        {
          numeric_type const str_weight(i_str.weight());
          if( str_weight == 0e0 ) continue; // don't fit this structure.

          t_FittingSet::str_iterator::rep_iterator i_rep = i_str.begin();
          t_FittingSet::str_iterator::rep_iterator const i_rep_end = i_str.end();
          t_Values::str_iterator::rep_iterator i_rep_val = i_str_val.begin();
          vector_type G(vec_size, 0); // values to sum to matrix and vector.

          // loop over structure representations.
          for(; i_rep != i_rep_end; ++i_rep, ++i_rep_val )
          {
            numeric_type const rep_weight(i_rep.weight() * str_weight);

            typedef t_FittingSet::str_iterator::rep_iterator
                                ::coordinate_iterator coordinate_iterator;
            typedef t_Values::str_iterator::rep_iterator
                            ::rank_iterator rank_iterator;
            rank_iterator i_rank_val = i_rep_val.begin();
            rank_iterator const i_rank_val_end = i_rep_val.end();
            t_ScalingFactors::const_iterator i_scale = scaling_factors_.begin();
#           ifdef LADA_DEBUG
              t_ScalingFactors::const_iterator const i_scale_end = scaling_factors_.end();
#           endif
            coordinate_iterator i_coord = i_rep.begin();


            // loop over ranks.
            for(size_t i(0); i_rank_val != i_rank_val_end; ++i_rank_val, ++i_scale )
            {
              LADA_ASSERT(i_scale != i_scale_end, "Iterator out of range.\n")
              numeric_type factor_i( i_rank_val.other() * (*i_scale) );
              rank_iterator::function_iterator i_func = i_rank_val.begin();
              rank_iterator::function_iterator const i_func_end = i_rank_val.end();
              // loop over inner functions.
              for(; i_func != i_func_end; ++i_func, i+=Functions::N, ++i_coord)
                G( i + (*i_coord) ) += rep_weight * (*i_func) * factor_i;
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
        for( size_t i(0); i < vec_size; ++i )
          for( size_t j(i+1); j < vec_size; ++j )
            _matrix(i,j) = _matrix(j,i);
      }

      void Collapse::add(Crystal::TStructure<std::string> const &_structure )
      {
        size_t Nmax(0);
        SumOfSeparables::const_iterator i_first = sumofseps_.begin();
        SumOfSeparables::const_iterator i_end = sumofseps_.end();
        for(; i_first != i_end; ++i_first)
          Nmax = std::max(Nmax, i_first->get<0>().size());
           
        LADA_ASSERT( (Nmax+2) % 3 == 0, "Unexpected number of variables.\n" )
        Representation representation(_structure, (Nmax + 2) / 3 );
        fitting_set_.add( representation, _structure.energy, _structure.weight );
        values_.add( representation, sumofseps_ );
      }
    } // namespace collapse
  } // namespace atomic_potential
} // namespace LaDa
