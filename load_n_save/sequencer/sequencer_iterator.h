//
//  Version: $Id: sequencer_iterator.h 1266 2009-08-10 05:01:26Z davezac $
//
#ifndef _LADA_XPR_SEQUENCE_ITERATOR_H_
#define _LADA_XPR_SEQUENCE_ITERATOR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/debug.h>
#include "sequencer.h"

namespace LaDa
{
  namespace load_n_save
  {
    namespace sequencer
    {
      //! Type of the iterator over ranges of and operation.
      template< class T_ITERATOR >
        struct const_iterator
        {
          public:
            //! Type of the values.
            typedef typename T_ITERATOR::value_type value_type;
            //! Type of the reference.
            typedef typename T_ITERATOR::reference reference;
            //! Type of the reference.
            typedef typename T_ITERATOR::pointer pointer;
        
            //! Constructor.
            const_iterator( Sequence::const_iterator const &_seq, T_ITERATOR const &_out )
                     : i_seq_(_seq), i_out_(_out) {}
            //! Constructor.
            const_iterator( const_iterator const &_c )
                     : i_seq_(_c.i_seq_), i_out_(_c.i_out_) {}
        
            //! Deref.
            reference operator*() const { return *i_out_; }
            //! Deref.
            pointer operator->() const { return &(*i_out_); }
            //! Pre-increment operator.
            const_iterator& operator++()
            {
              if( is_object() ) ++i_out_;
              ++i_seq_;
              return *this;
            }
            //! Post-increment operator.
            const_iterator operator++(int) { return const_iterator( ++(*this) ); }
            //! Returns true if iterators are at same position.
            bool operator==( const_iterator const &_c ) const
            {
#             ifdef _LADADEBUG
                if( i_seq_ == _c.i_seq_  ) 
                  { LADA_ASSERT( i_out_ == _c.i_out_, "Incoherent iterators.\n" ); }
#             endif
              return i_seq_ == _c.i_seq_;
            }
            //! Returns true if iterators are not at same position.
            bool operator!=( const_iterator const &_c ) const { return not( *this == _c ); }
            //! True if tag is an "or" operator.
            bool is_or() const { return *i_seq_ == or_; }
            //! True if tag is an object.
            bool is_object() const { return *i_seq_ == object; }
            //! True if tag starts a group.
            bool is_start_group() const { return *i_seq_ == start_group; }
            //! True if tag ends a group.
            bool is_end_group() const { return *i_seq_ == end_group; }
          protected:
            //! Holds an iterator to the sequence definition.
            Sequence::const_iterator i_seq_;
            //! Holds an iterator to container of objects.
            T_ITERATOR i_out_;
        };

      namespace details
      {
        template< class T_ITERATOR >
          void find_group_end_impl( T_ITERATOR &_first, T_ITERATOR const &_end )
          {
            for(; _first != _end; ++_first )
              if( _first.is_start_group() )
              {
                ++_first;
                find_group_end_impl( _first, _end );
                LADA_DOASSERT( _first != _end, "Should not be at container end.\n" );
              }
              else if( _first.is_end_group() ) break;
          };
      }
      template< class T_ITERATOR > 
        void find_group_end( T_ITERATOR &_first, T_ITERATOR const &_end )
        {
          if( _first == _end ) return;
          if( not _first.is_start_group() ) return;
          details::find_group_end_impl( ++_first, _end );
        }

    }  // end of sequencer namespace
  } // end of load_n_save namespace.
} // namespace LaDa

#endif

