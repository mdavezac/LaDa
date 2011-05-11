//
//  Version: $Id: sequencer.h 1266 2009-08-10 05:01:26Z davezac $
//

#ifndef _LADA_LNS_XPR_SEQUENCER_H_
#define _LADA_LNS_XPR_SEQUENCER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <list>
#ifdef _LADADEBUG
# include<iostream>
#endif 

#include <opt/debug.h>

namespace LaDa 
{
  namespace load_n_save
  {
    //! \brief Helps transform a linear sequence into a sequence with operations.
    //! \details a linear sequence of objects a, b, c, d, ... can be
    //!          transformed into an algebraic expression (a and b) or c and d... using
    //!          this object. This object stores the sequence of & and | and
    //!          grouping (,) operators. The sequence of object must be held
    //!          elsewhere.
    namespace sequencer
    {
      //! Operator tags.
      enum tags
      {
        object      = 0, //!< just an object.
        start_group = 1, //!< starts a group with (.
        end_group   = 2, //!< ends a groupt with ).
        or_         = 3  //!< explicit | operator. & operators are implied. 
      };

      //! Type which holds the operator sequence.
      class Sequence : public std::list<tags>
      {
        public:
          Sequence() {};
          Sequence( Sequence &_sequence) : std::list<tags>(_sequence) {};
          Sequence( size_t _n ) : std::list<tags>(_n, object) {};
          Sequence( std::string const &_str ) { transform_(_str); }
          Sequence( char const *_str ) { transform_( std::string(_str) ); }

        private:
          void transform_( std::string const &_str );
      };

#     ifdef _LADADEBUG
        std::ostream& operator<<( std::ostream &_stream, Sequence const& _s );
#      define LADA_PRINT(a) std::cout << a << "\n";
#     else
#      define LADA_PRINT(a)
#     endif



      inline void operator|=( Sequence &_a, Sequence const &_b ) 
      {
        _a.push_back(or_); 
        std::copy( _b.begin(), _b.end(), std::back_inserter(_a) );
        LADA_PRINT(_a); 
      }
      void operator&=( Sequence &_a, Sequence const &_b );
#    undef LADA_PRINT
  
      //! \brief If *\a_first == start_group, returns iterator at close of parenthesis,
      //!       returns \a _first otherwise.
      template< class T_IT1, class T_IT2 >
        void find_group_end( T_IT1 &_first, T_IT1 const &_end, T_IT2 &_out );
  
  
      namespace details
      {
        template< class T_IT, class T_OUT >
          void find_group_end_impl( T_IT &_first, T_IT const &_end, T_OUT &_out )
          {
            for(; _first != _end; ++_first )
              if( *_first == start_group ) 
              {
                ++_first;
                find_group_end_impl( _first, _end, _out );
                LADA_DOASSERT( _first != _end, "Should not be at container end.\n" );
              }
              else if( *_first == end_group ) break;
              else if( *_first == object  ) ++_out;
          };
      }
  
      template< class T_IT, class T_OUT >
        void find_group_end( T_IT &_first, T_IT const &_end, T_OUT &_out )
        {
          if( _first == _end ) return;
          if( *_first != start_group ) return;
          details::find_group_end_impl( ++_first, _end, _out );
        }
      template< class T_IT, class T_OUT >
        void find_group_end( std::pair<T_IT,T_IT> &_seq, std::pair<T_OUT, T_OUT> &_cont )
          { find_group_end( _seq.first, _seq.second, _cont.first  ); }

      template< class T_IT >
        bool is_next_or( T_IT const &_first, T_IT const &_end )
        {
          if( _first == _end ) return false; 
          T_IT first( _first + 1 );
          if( first == _end ) return false; 
          return *first == or_;
        }

      bool is_all_ands( Sequence const &_a );

    } // namespace sequencer.
  } // namespace load_n_save
} // namespace LaDa

#endif
