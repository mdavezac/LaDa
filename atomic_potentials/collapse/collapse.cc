//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/vector_proxy.hpp>

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
            for(size_t i(0); rank_range; ++rank_range)
            {
              const_iterator i_first = rank_range.begin();
              const_iterator const i_end = rank_range.end();
              for(; i_first != i_end; ++i_first)
                for( size_t j(0); j < Functions::N; ++j, ++i)
                  vec(i) = (*i_first)[j];
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
      };

      numeric_type Collapse::convergence() const
      {
        numeric_type result(0);
        size_t nstr(0);
        t_FittingSet :: str_iterator i_str = fitting_set_.begin(0);
        t_FittingSet :: str_iterator const i_str_end = fitting_set_.end(0);
        t_Values::const_str_iterator i_str_val = values_.begin(0);
        //! Loop over structures.
        for(; i_str != i_str_end; ++i_str, ++i_str_val)
        {
          numeric_type const str_weight(i_str.weight());
          if( str_weight == 0e0 ) continue; // don't fit this structure.
          ++nstr;
        
          t_FittingSet::str_iterator::rep_iterator i_rep = i_str.begin();
          t_FittingSet::str_iterator::rep_iterator const i_rep_end = i_str.end();
          t_Values::const_str_iterator::const_rep_iterator i_rep_val = i_str_val.begin();
        
          // loop over structure representations.
          numeric_type str_val(-i_str.energy()); 
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
          result += str_weight * str_val;
        } // loop over structures.
        return result / numeric_type(nstr); 
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

      size_t Collapse::nb_coordinates() const { return sumofseps_.nb_coordinates(); }



    } // namespace collapse
  } // namespace atomic_potential
} // namespace LaDa
