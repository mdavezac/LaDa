//
//  Version: $Id: binary_range.h 1266 2009-08-10 05:01:26Z davezac $
//
#ifndef _LADA_LNS_SEQUENCER_BINARY_RANGE_H_
#define _LADA_LNS_SEQUENCER_BINARY_RANGE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/debug.h>
#include "binary.h"
#include "binary_const_iterator.h"

namespace LaDa
{
  namespace load_n_save
  {
    namespace sequencer
    {
      //! Type of the iterator over ranges of and operation.
      template< class T_FIRST, class T_SECOND >
        struct BinaryRange : public std::pair
          < 
            binary_const_iterator<T_FIRST, T_SECOND>,
            binary_const_iterator<T_FIRST, T_SECOND> 
          >
        {
          public:
            //! Type of the internal iterators.
            typedef binary_const_iterator<T_FIRST, T_SECOND> value_type;
            //! Type of the base class
            typedef std::pair<value_type, value_type> t_Base;
            //! Constructor.
            template< class T_FCONTAINER, class T_SCONTAINER >
              BinaryRange   ( Binary& _seq, T_FCONTAINER &_first, T_SCONTAINER &_second )
                          : t_Base( value_type(_seq.begin(), _first.begin(), _second.begin()),
                                    value_type(_seq.begin(), _first.begin(), _second.begin()) ),
                            end_(_seq.end(), _first.end(), _second.end()) { ++(*this); }
            //! Constructor.
            BinaryRange   ( Binary::const_iterator const &_seq0, Binary::const_iterator const &_seq1, 
                            T_FIRST const &_first, T_FIRST const &_first_end,
                            T_SECOND const &_second, T_SECOND const &_second_end )
                       : t_Base( value_type(_seq0, _first, _second), value_type(_seq0, _first, _second) ),
                         end_(_seq1, _first_end, _second_end) { ++(*this); }
            //! Constrctor.
            BinaryRange   ( Binary const &_seq, 
                            T_FIRST const &_first, T_FIRST const &_first_end,
                            T_SECOND const &_second, T_SECOND const &_second_end )
                        : t_Base( value_type(_seq.begin(), _first, _second),
                                  value_type(_seq.begin(), _first, _second) ),
                          end_(_seq.end(), _first_end, _second_end)  { ++(*this); }
            //! Segment Constructor.
            BinaryRange   ( value_type &_first, value_type &_end )
                        : t_Base(_first, _first), end_(_end) { ++(*this); }
            //! Copy Constructor.
            BinaryRange   ( BinaryRange const &_c )
                        : t_Base(_c), end_(_c.end_) {}
        
            //! Pre-increment operator.
            BinaryRange& operator++();
            //! Post-increment operator.
            BinaryRange operator++(int) { return BinaryRange( ++(*this) ); }
            //! Returns true if not end of sequence.
            operator bool() const { return first != end_; }

            using t_Base::first;
            using t_Base::second;
          protected:
            //! End of the sequence.
            value_type end_;
        };

      //! Type of this iterator for constant containers.
      template< class T_FCONT, class T_SCONT >
        struct binary_const_range
        {
          //! Iteratortype.
          typedef BinaryRange<typename T_FCONT::const_iterator, typename T_SCONT::const_iterator> type;
        }; 

      //! Type of the iterator over ranges of and operation.
      template< class T_FIRST, class T_SECOND >
        BinaryRange<T_FIRST, T_SECOND>& BinaryRange<T_FIRST, T_SECOND> :: operator++()
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

