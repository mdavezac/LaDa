//
//  Version: $Id: binary_const_iterator.h 1266 2009-08-10 05:01:26Z davezac $
//
#ifndef _LADA_LNS_SEQUENCER_BINARY_CONST_ITERATOR_H_
#define _LADA_LNS_SEQUENCER_BINARY_CONST_ITERATOR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/debug.h>
#include "binary.h"

namespace LaDa
{
  namespace load_n_save
  {
    namespace sequencer
    {
      //! Type of the iterator over ranges of and operation.
      template< class T_FIRST, class T_SECOND >
        struct binary_const_iterator
        {
          public:
            //! Type of the values of the first iterator.
            typedef typename T_FIRST::value_type first_value_type;
            //! Type of a reference to the values of the first iterator.
            typedef typename T_FIRST::reference first_reference;
            //! Type of the pointer to the values of the first iterator.
            typedef typename T_FIRST::pointer first_pointer;
            //! Type of the values of the second iterator.
            typedef typename T_SECOND::value_type second_value_type;
            //! Type of a reference to the values of the second iterator.
            typedef typename T_SECOND::reference second_reference;
            //! Type of the pointer to the values of the second iterator.
            typedef typename T_SECOND::pointer second_pointer;
        
            //! Constructor.
            binary_const_iterator   ( Binary::const_iterator const &_seq, 
                                      T_FIRST const &_first,
                                      T_SECOND const &_second )
                                  : i_seq_(_seq), i_first_(_first), i_second_(_second) {}
            //! Constructor.
            binary_const_iterator   ( binary_const_iterator const &_c )
                                  : i_seq_(_c.i_seq_), i_first_(_c.i_first_), i_second_(_c.i_second_) {}
        
            //! Pointer to first values.
            first_pointer first() const { return &(*i_first_); }
            //! Pointer to second values.
            second_pointer second() const { return &(*i_second_); }
            //! Pre-increment operator.
            binary_const_iterator& operator++()
            {
              if( is_first() ) ++i_first_;
              else if( is_second() ) ++i_second_;
              ++i_seq_;
              return *this;
            }
            //! Post-increment operator.
            binary_const_iterator operator++(int) { return binary_const_iterator( ++(*this) ); }
            //! Returns true if iterators are at same position.
            bool operator==( binary_const_iterator const &_c ) const
            {
#             ifdef _LADADEBUG
                if( i_seq_ == _c.i_seq_  )
                {
                  LADA_ASSERT( i_first_  ==  _c.i_first_, "Incoherent iterators.\n" ); 
                  LADA_ASSERT( i_second_ == _c.i_second_, "Incoherent iterators.\n" ); 
                }
#             endif
              return i_seq_ == _c.i_seq_;
            }
            //! Returns true if iterators are not at same position.
            bool operator!=( binary_const_iterator const &_c ) const { return not( *this == _c ); }
            //! True if tag is an "or" operator.
            bool is_or() const { return *i_seq_ == or_; }
            //! True if tag is a first object.
            bool is_first() const { return *i_seq_ == first_object; }
            //! True if tag is a second object.
            bool is_second() const { return *i_seq_ == second_object; }
            //! True if tag starts a group.
            bool is_start_group() const { return *i_seq_ == start_group; }
            //! True if tag ends a group.
            bool is_end_group() const { return *i_seq_ == end_group; }
          protected:
            //! Holds an iterator to the sequence definition.
            Binary::const_iterator i_seq_;
            //! Holds an iterator to container of first objects.
            T_FIRST i_first_;
            //! Holds an iterator to container of second objects.
            T_SECOND i_second_;
        };

      namespace details
      {
        template< class T_FIRST, class T_SECOND >
          void find_group_end_impl( binary_const_iterator<T_FIRST, T_SECOND> &_first,
                                    binary_const_iterator<T_FIRST, T_SECOND> const &_end )
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
      template< class T_FIRST, class T_SECOND >
        void find_group_end( binary_const_iterator<T_FIRST, T_SECOND> &_first,
                             binary_const_iterator<T_FIRST, T_SECOND> const &_end )
        {
          if( _first == _end ) return;
          if( not _first.is_start_group() ) return;
          details::find_group_end_impl( ++_first, _end );
        }

    }  // end of sequencer namespace
  } // end of load_n_save namespace.
} // namespace LaDa

#endif

