//
//  Version: $Id$
//
#ifndef LADA_ATOMIC_POTENTIAL_VARIABLE_MAJOR_H_
#define LADA_ATOMIC_POTENTIAL_VARIABLE_MAJOR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/fusion/include/at_c.hpp>

#include "../sum_of_separables.h"

namespace LaDa
{
  namespace atomic_potential
  {
    namespace collapse
    {
      //! Inverts iteration, such that it goes over variables first, and ranks second.
      class CoordRange
      {
#       ifdef LADA_VMRASSERT
#         error LADA_VMRASSERT already defined.
#       endif
#       ifdef LADA_DEBUG
#         define LADA_VMRASSERT LADA_ASSERT( index_ < max_index_, "Range out of range.\n") 
#       else
#         define LADA_VMRASSERT 
#       endif
        public:
          //! Type of the dereferenced value.
          typedef size_t value_type;

          //! Iterator over ranks.
          class rank_range;
          //! Iterator over ranks.
          class const_rank_range;
          //! Constructor.
          CoordRange   (SumOfSeparables const &_sumofseps )
                     : sumofseps_(_sumofseps), index_(0),
                       max_index_(sumofseps_.nb_coordinates() ) {}
          //! Copy Constructor.
          CoordRange   (CoordRange const &_c)
                     : sumofseps_(_c.sumofseps_),
                       index_(_c.index_), max_index_(_c.max_index_) {}
      
          //! Range over ranks.
          rank_range range();
          //! Range over ranks.
          const_rank_range range() const;

          //! Increments.
          bool operator++() { ++index_; return this->operator bool(); }

          //! Continue iterating or not.
          operator bool() const { return index_ >= max_index_; }

          //! Returns size of the range.
          size_t size() const { return max_index_; }

          value_type operator*() const { return index_; }

        private:
          //! Type of the numeric values.
          SumOfSeparables const &sumofseps_;
          //! Current index.
          size_t index_;
          //! Maximum index.
          size_t max_index_;

      };

#     ifdef LADA_RANKRANGE
#       error LADA_RANKRANGE already defined.
#     endif
#     define LADA_RANKRANGE 0
#     include "variable_major.rank_range.h"

      inline CoordRange::rank_range CoordRange::range() 
      {
        LADA_VMRASSERT; 
        return rank_range( sumofseps_.begin(), sumofseps_.end(), index_ ); 
      }
      inline CoordRange::const_rank_range CoordRange::range() const
      {
        LADA_VMRASSERT; 
        return const_rank_range( sumofseps_.begin(), sumofseps_.end(), index_ ); 
      }
#     undef LADA_VMRASSERT


    } // namspace collapse
  } // namespace atomic_potential
} // namespace LaDa
#endif
