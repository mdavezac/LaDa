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
          //! Constructor.
          CoordRange(SumOfSeparables const &_sumofseps );
          //! Copy Constructor.
          CoordRange   (CoordRange const &_c)
                     : sumofseps_(_c.sumofseps_),
                       index_(_c.index_), max_index_(_c.max_index_) {}
      
          //! Iterator over ranks.
          rank_range begin() const;
          //! Iterator over ranks.
          rank_range end() const;

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

      //! Compares two rank iterators.
      bool operator==( CoordRange::rank_range const& _a,
                       CoordRange::rank_range const& _b);

      class CoordRange :: rank_range
      {
        friend bool operator==( rank_range const &, rank_range const& );
        friend class CoordRange;
        //! Type of the var iterator.
        typedef SumOfSeparables::t_Function::const_iterator var_iterator_;
        public:
          //! value type
          typedef SumOfSeparables::t_Function::const_iterator::value_type value_type;
          //! reference type
          typedef SumOfSeparables::t_Function::const_iterator::reference reference;
          //! pointer type
          typedef SumOfSeparables::t_Function::const_iterator::pointer pointer;
          
          //! Iterator over functions and coefficients for set rank and variable.
          typedef SumOfSeparables::t_Function::t_Function::const_iterator inner_iterator;

          //! Copy Constructor.
          rank_range   (rank_range const &_c)
                        : i_rank_(_c.i_rank_), i_rank_end_(_c.i_rank_end_), index_(_c.index_) {}
      
          //! Iterator over ranks.
          inner_iterator begin() const { return get_i_var_()->begin(); }
          //! Iterator over ranks.
          inner_iterator end() const { return get_i_var_()->end(); }

          //! Increments.
          void operator++();

          //! Deref.
          value_type operator*() const { return get_i_var_().operator*(); }
          //! Deref.
          pointer operator->() const { return get_i_var_().operator->(); }
          //! True if still iteratable.
          operator bool() const { return i_rank_ != i_rank_end_; }

        private:
          //! Constructor.
          rank_range   (SumOfSeparables::const_iterator const &_irank,
                        SumOfSeparables::const_iterator const &_irank_end,
                        size_t _i);

          //! Gets correct variable iterator.
          var_iterator_ get_i_var_() const;

          //! Current rank positions.
          SumOfSeparables::const_iterator i_rank_;
          //! Last rank position.
          SumOfSeparables::const_iterator i_rank_end_;
          //! Variable index.
          size_t index_;
      };

      inline CoordRange::rank_range CoordRange::begin() const 
        { LADA_VMRASSERT; return rank_range( sumofseps_.begin(), sumofseps_.end(), index_ ); }
          //! Iterator over ranks.
      inline CoordRange::rank_range CoordRange::end() const
        { LADA_VMRASSERT; return rank_range( sumofseps_.end(), sumofseps_.end(), index_ ); }
#       undef LADA_VMRASSERT

      inline bool operator==( CoordRange::rank_range const& _a,
                              CoordRange::rank_range const& _b)
      {
        LADA_ASSERT( _a.index_ == _b.index_, "Inconsistent iterators.\n")
        return _a.i_rank_ == _b.i_rank_; 
      }
      //! Compares two rank iterators.
      inline bool operator!=( CoordRange::rank_range const& _a,
                              CoordRange::rank_range const& _b)
        { return not (_a == _b); }

    } // namspace collapse
  } // namespace atomic_potential
} // namespace LaDa
#endif
