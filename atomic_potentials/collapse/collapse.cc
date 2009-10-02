//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cmath>
#include <numeric>
#include <algorithm>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/lambda/lambda.hpp>

#include <crystal/structure.h>

#include "../sum_of_separables.h"
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
      Collapse::Collapse   (SumOfSeparables & _sumofseps) 
                         : sumofseps_(_sumofseps)
      {
        try
        {
          typedef SumOfSeparables::const_coord_range const_coord_range;
          typedef const_coord_range::const_rank_range const_rank_range;
          typedef const_rank_range::const_iterator const_iterator;
          coefficients_.clear();
          nb_funcs_.clear(); nb_funcs_.reserve( _sumofseps.nb_coordinates() );
          for(const_coord_range range( _sumofseps.const_range() ); range; ++range)
          {
            size_t nmax(0);
            const_coord_range::const_rank_range rank_range(range.range());
            for(; rank_range; ++rank_range)
            {
              const_iterator i_first = rank_range.begin();
              const_iterator const i_end = rank_range.end();
              for(; i_first != i_end; ++i_first, nmax += Functions::N);
            }
          
            vector_type vec(nmax);
            rank_range = range.range();
            nb_funcs_.resize(nb_funcs_.size() + 1);
            for(size_t i(0); rank_range; ++rank_range)
            {
              const_iterator i_first = rank_range.begin();
              const_iterator const i_end = rank_range.end();
              size_t n(0);
              for(; i_first != i_end; ++i_first, ++n)
                for( size_t j(0); j < Functions::N; ++j, ++i)
                  vec(i) = (*i_first)[j];
              nb_funcs_.back().push_back(n);
            }
            coefficients_.push_back(vec);
          }

          SumOfSeparables::const_iterator i_scale( _sumofseps.begin() );
          SumOfSeparables::const_iterator const i_scale_end( _sumofseps.end() );
          for(; i_scale != i_scale_end; ++i_scale)
            scaling_factors_.push_back( i_scale->get<1>() );

        }
        catch(...)
        {
          std::cerr << "Could not create collapse object.\n";
        }
      }

      void Collapse::reassign() const
      {
        // coefficients
        typedef SumOfSeparables::coord_range t_coord_range;
        typedef t_coord_range::rank_range t_rank_range;
        typedef t_rank_range::iterator iterator;
        for(t_coord_range range( sumofseps_.range() ); range; ++range)
        {
          vector_type const &vec(coefficients_[*range]);
          t_rank_range rank_range = range.range();
          for(size_t i(0); rank_range; ++rank_range)
          {
            iterator i_first = rank_range.begin();
            iterator const i_end = rank_range.end();
            for(; i_first != i_end; ++i_first)
              for( size_t j(0); j < Functions::N; ++j, ++i)
                (*i_first)[j] = vec(i);
          }
        }

        // scales
        LADA_ASSERT( sumofseps_.size() == scaling_factors_.size(), "Incoherent containers." )
        SumOfSeparables::iterator i_rank = sumofseps_.begin();
        SumOfSeparables::iterator const i_rank_end = sumofseps_.end();
        t_ScalingFactors::const_iterator i_scale = scaling_factors_.begin();
        for(; i_rank != i_rank_end; ++i_rank, ++i_scale)
          boost::tuples::get<1>(*i_rank) = *i_scale;
      };

      // Removes a header dependence by defining it here.
      size_t Collapse :: nb_coordinates() const { return sumofseps_.nb_coordinates(); }

      numeric_type Collapse::convergence() const
      {
        numeric_type result(0);
        numeric_type nstr(0);
        t_FittingSet :: str_iterator i_str = fitting_set_.begin(0);
        t_FittingSet :: str_iterator const i_str_end = fitting_set_.end(0);
        t_Values::const_str_iterator i_str_val = values_.begin(0);
        //! Loop over structures.
        for(; i_str != i_str_end; ++i_str, ++i_str_val)
        {
          numeric_type const str_weight(i_str.weight());
          if( str_weight == 0e0 ) continue; // don't fit this structure.
          nstr += str_weight;
        
          t_FittingSet::str_iterator::rep_iterator i_rep = i_str.begin();
          t_FittingSet::str_iterator::rep_iterator const i_rep_end = i_str.end();
          t_Values::const_str_iterator::const_rep_iterator i_rep_val = i_str_val.begin();
        
          // loop over structure representations.
          numeric_type str_val(0);
          for(; i_rep != i_rep_end; ++i_rep, ++i_rep_val )
          {
            numeric_type const rep_weight(i_rep.weight() * str_weight);
        
            typedef t_Values::const_str_iterator::const_rep_iterator
                            ::const_rank_iterator rank_iterator;
            rank_iterator i_rank_val = i_rep_val.begin();
            rank_iterator const i_rank_val_end = i_rep_val.end();
            t_ScalingFactors::const_iterator i_scale = scaling_factors_.begin();
#           ifdef LADA_DEBUG
              t_ScalingFactors::const_iterator const i_scale_end = scaling_factors_.end();
#           endif
        
            // loop over ranks.
            for(size_t i(0); i_rank_val != i_rank_val_end; ++i_rank_val, ++i_scale )
            {
              LADA_ASSERT(i_scale != i_scale_end, "Iterator out of range.\n");
              str_val += i_rank_val.all() * (*i_scale) * rep_weight;
            } // loop over ranks
          } // loop over representations
          std::cout << "str_val: " << str_val << "\n";
          result += str_weight * (str_val - i_str.energy()) * (str_val - i_str.energy());
        } // loop over structures.
        return result / nstr; 
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

      void Collapse::update(size_t _i, vector_type const &_x)
      { 
        namespace bl = boost::lambda;
        LADA_DOASSERT( _i < nb_funcs_.size(), "Index out-of-range.\n" )
        LADA_ASSERT( coefficients_.size() == nb_funcs_.size(), "Incoherent containers.\n" )
        LADA_ASSERT( scaling_factors_.size() == nb_funcs_[_i].size(), "Incoherent containers.\n" )

        coefficients_[_i] = _x;
        values_.update(coefficients_[_i], fitting_set_, _i);

        // Normalizes functionals.
        t_ScalingFactors :: iterator i_scale = scaling_factors_.begin();
        std::vector<size_t>::const_iterator i_nbfuncs = nb_funcs_[_i].begin();
        std::vector<size_t>::const_iterator const i_nbfuncs_end = nb_funcs_[_i].end();
        t_Coefficients::value_type::iterator i_coef = coefficients_[_i].begin();
#       ifdef LADA_DEBUG
          size_t acc(0);
#       endif
        for(; i_nbfuncs != i_nbfuncs_end; ++i_nbfuncs, ++i_scale) 
        {
          t_Coefficients::value_type::iterator const i_first = i_coef;
          i_coef += Functions::N * (*i_nbfuncs);
#         ifdef LADA_DEBUG
            acc += Functions::N * (*i_nbfuncs);
            LADA_ASSERT(acc <= coefficients_[_i].size(), "index out of range.\n") 
#         endif
          numeric_type const norm
          (
            std::sqrt( std::accumulate(i_first, i_coef, numeric_type(0), bl::_1 + bl::_2*bl::_2) )
          );
          numeric_type const inv_norm( numeric_type(1) / norm );
          std::for_each( i_first, i_coef, bl::_1 *= bl::constant(inv_norm) );
          (*i_scale) *= norm;
        }
      }

      numeric_type Collapse::y_squared() const
      {
        numeric_type result(0);
        t_FittingSet :: str_iterator i_str = fitting_set_.begin(0);
        t_FittingSet :: str_iterator const i_str_end = fitting_set_.end(0);
        //! Loop over structures.
        for(; i_str != i_str_end; ++i_str) result += i_str.weight() * i_str.energy();
        return result;
      }
      numeric_type Collapse::sum_w() const
      {
        numeric_type result(0);
        t_FittingSet :: str_iterator i_str = fitting_set_.begin(0);
        t_FittingSet :: str_iterator const i_str_end = fitting_set_.end(0);
        //! Loop over structures.
        for(; i_str != i_str_end; ++i_str) result += i_str.weight();
        return result;
      }

    } // namespace collapse
  } // namespace atomic_potential
} // namespace LaDa
