//
//  Version: $Id: sequencer_range.h 1266 2009-08-10 05:01:26Z davezac $
//
#ifndef _LADA_XPR_SEQUENCE_RANGE_ITERATOR_H_
#define _LADA_XPR_SEQUENCE_RANGE_ITERATOR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/debug.h>
#include "sequencer.h"
#include "sequencer_iterator.h"

namespace LaDa
{
  namespace load_n_save
  {
    namespace sequencer
    {
      //! Type of the iterator over ranges of and operation.
      template< class T_ITERATOR >
        struct AndRange : public std::pair< const_iterator<T_ITERATOR>, const_iterator<T_ITERATOR> >
        {
          public:
            //! Type of the internal iterators.
            typedef const_iterator<T_ITERATOR> value_type;
            //! Type of the base class
            typedef std::pair<value_type, value_type> t_Base;
            //! Constructor.
            template< class T_CONTAINER >
              AndRange   ( Sequence& _seq, T_CONTAINER &_cont )
                       : t_Base( value_type(_seq.begin(), _cont.begin()),
                                 value_type(_seq.begin(), _cont.begin()) ),
                         end_(_seq.end(), _cont.end()) { ++(*this); }
            //! Constructor.
            AndRange   ( Sequence::const_iterator const &_seq0, Sequence::const_iterator const &_seq1, 
                         T_ITERATOR const &_it0, T_ITERATOR const &_it1 )
                     : t_Base( value_type(_seq0, _it0), value_type(_seq0, _it0) ),
                       end_(_seq1, _it1) { ++(*this); }
            //! Constrctor.
            AndRange   ( Sequence const &_seq, T_ITERATOR const &_it0, T_ITERATOR const &_it1 )
                     : t_Base( value_type(_seq.begin(), _it0), value_type(_seq.begin(), _it0) ),
                       end_(_seq.end(), _it1)  { ++(*this); }
            //! Segment Constructor.
            AndRange   ( value_type &_first, value_type &_end )
                     : t_Base(_first, _first), end_(_end) { ++(*this); }
            //! Copy Constructor.
            AndRange   ( AndRange const &_c )
                     : t_Base(_c), end_(_c.end_) {}
        
            //! Pre-increment operator.
            AndRange& operator++();
            //! Post-increment operator.
            AndRange operator++(int) { return AndRange( ++(*this) ); }
            //! Returns true if not end of sequence.
            operator bool() const { return first != end_; }

            using t_Base::first;
            using t_Base::second;
          protected:
            //! End of the sequence.
            value_type end_;
        };

      //! Type of this iterator for a constant container.
      template< class T_CONTAINER >
        struct const_range_iterator
        {
          //! Iteratortype.
          typedef AndRange<typename T_CONTAINER::const_iterator> type;
        }; 

      //! Type of the iterator over ranges of and operation.
      template< class T_ITERATOR >
        AndRange<T_ITERATOR>& AndRange<T_ITERATOR> :: operator++()
        { 
          first = second;
          if( first.is_or() ) { ++first; ++second; }
 
          for(; second != end_; ++second )
          {
            if( second.is_or() ) break;
            else if( second.is_start_group() ) find_group_end( second, end_ );
#           ifdef _LADADEBUG
            else if( second.is_end_group() ) { __THROW_ERROR("Bug -- Should not reach this point.\n"); }
#           endif
          }
 
          return *this;
        }

    }  // end of sequencer namespace
  } // end of load_n_save namespace.
} // namespace LaDa

#endif

