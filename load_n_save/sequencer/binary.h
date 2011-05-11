
//
//  Version: $Id: binary.h 1266 2009-08-10 05:01:26Z davezac $
//

#ifndef _LADA_LNS_SEQUENCER_BINARY_H_
#define _LADA_LNS_SEQUENCER_BINARY_H_

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
        first_object  = 0, //!< just an object.
        second_object = 1, //!< just another object.
        start_group   = 2, //!< starts a group with (.
        end_group     = 3, //!< ends a groupt with ).
        or_           = 4  //!< explicit | operator. & operators are implied. 
      };

      //! Helps convert two seqences into binary trees.
      class Binary : public std::list<tags>
      {
        public:
          //! Constructor.
          Binary() {};
          //! Copy Constructor.
          Binary( Binary &_sequence) : std::list<tags>(_sequence) {};
          //! Constructor.
          Binary( std::string const &_str ) { transform_(_str); }
          //! Constructor.
          Binary( char const *_str ) { transform_( std::string(_str) ); }

        private:
          //! Transforms a string into a sequence.
          void transform_( std::string const &_str );
      };

#     ifdef _LADADEBUG
        std::ostream& operator<<( std::ostream &_stream, Binary const& _s );
#      define LADA_PRINT(a) std::cout << a << "\n";
#     else
#      define LADA_PRINT(a)
#     endif



      //! Joins two sequences via an || operator.
      inline void operator|=( Binary &_a, Binary const &_b ) 
      {
        _a.push_back(or_); 
        std::copy( _b.begin(), _b.end(), std::back_inserter(_a) );
        LADA_PRINT(_a); 
      }
      //! Joins two sequences via an && operator.
      void operator&=( Binary &_a, Binary const &_b );
#    undef LADA_PRINT
  
      //! \brief If *\a_first == start_group, returns iterator at close of parenthesis,
      //!       returns \a _first otherwise.
      template< class T_IT1, class T_IT2, class T_IT3 >
        void find_group_end( T_IT1 &_first, T_IT1 const &_end, T_IT2 &_out, T_IT3 &_out2 );
  
  
      namespace details
      {
        template< class T_IT, class T_OUT, class T_OUT2 >
          void find_group_end_impl( T_IT &_first, T_IT const &_end, T_OUT &_out, T_OUT2 &_out2 )
          {
            for(; _first != _end; ++_first )
              if( *_first == start_group ) 
              {
                ++_first;
                find_group_end_impl( _first, _end, _out, _out2 );
                LADA_DOASSERT( _first != _end, "Should not be at container end.\n" );
              }
              else if( *_first == end_group ) break;
              else if( *_first == first_object  ) ++_out;
              else if( *_first == second_object  ) ++_out2;
          };
      }
  
      template< class T_IT, class T_OUT, class T_OUT2 >
        void find_group_end( T_IT &_first, T_IT const &_end, T_OUT &_out, T_OUT2 &_out2 )
        {
          if( _first == _end ) return;
          if( *_first != start_group ) return;
          details::find_group_end_impl( ++_first, _end, _out );
        }

      bool is_all_ands( Binary const &_a );

    } // namespace sequencer.
  } // namespace load_n_save
} // namespace LaDa

#endif
