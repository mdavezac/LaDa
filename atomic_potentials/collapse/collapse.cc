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
