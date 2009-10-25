//
//  Version: $Id$
//
#ifndef LADA_OPT_NDIMITERATOR_H_
#define LADA_OPT_NDIMITERATOR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <vector>
#include <algorithm>

#include "opt/types.h"


namespace LaDa
{
  namespace opt
  {

    /** \brief Iterator of dynamic dimensionality for loop of dynamic nestedness.
     * \details Creates an iterator which allows for a dynamic number of loops whithin loops. 
     *          Imagine for instance that we thant a three-times nested loop
     *          where each loops iterates 4 times, then we could Ndim_Iterator in
     *          the following code.
       \code 
         // Creates a dimensional iterator of integers with operator <= to check
         // for the conditional end of a loop.
         opt::Ndim_Iterator< types::t_int, std::less_equal<types::t_int> > iterator;
         // adds outermost loop
         iterator.add( 0, 3);
         // adds intermediate loop
         iterator.add( 0, 3);
         // adds innermost loop
         iterator.add( 0, 3);

         do
         {
           // value of the outer loop
           types::t_int outer        =  iterator.access(0);
           // value of the intermediate loop
           types::t_int intermediate =  iterator.access(1);
           // value of the innner loop
           types::t_int inner        =  iterator.access(2);

           ...


         } while( ++iterator );
        \endcode
                This code is equivalent to the following
        \code
          for( types::t_int outer = 0; outer <= 3; ++outer )
            for( types::t_int intermediate = 0; intermediate <= 3; ++intermediate )
              for( types::t_int inner = 0; inner <= 3; ++innner )
              {
                 ...
              }
        \endcode
                The advantage of the second code is that it is more compact ;).
                The advantage of the first code is that we do not need to know in
                advance how nested the loop will be. Indeed, rather than
                explicitely adding three dimensions to the N-dimensional
                iterators, we could add as... whatever crops up. 

                In any case, the N-dimensional iterators needs all the following
                to construct a loop 
        \code 
           for( T_OBJECT object(_start); _condition( object, _end ); ++object )
        \endcode
                It may not be obvious, but here is the input to the loop:
                  - T_OBJECT: is a type defining a quantity
                  - _condition: is a binary functor which returns true as long as
                                the loop should continue.
                  _ operator++(): T_OBJECT understands the incrementation operator.
                  _ operator=( ? ): T_OBJECT must have a constructor whcih can
                                    take _start as an argument.
                  - _start: an initial value for the object
                  - _end: an end value for the object
                  .
                The first two items are the template arguments to Ndim_Iterator.
                The second two are behaviors required by Ndim_Iterator for
                T_OBJECT. Note that _condition may itself require some behaviors
                of T_OBJECT. The last two items are simply meant to initialize
                and end the loop.

                Finally, other than the member function Ndim_Iterator::access()
                displayed above, one can also a dynamic loop over inner
                iterators, where
         \code
           // value of the outer loop
           types::t_int outer        =  iterator.access(0);
           // value of the intermediate loop
           types::t_int intermediate =  iterator.access(1);
           // value of the innner loop
           types::t_int inner        =  iterator.access(2);
         \endcode
                can become,
         \code
           std::vector< T_OBJECT > inners;
           for( iterator.init_loop(); iterator.loop_is_valid(); iterator.next_iterator() );
             inners.push_back( iterator.get_current_iterator() );
         \endcode
                where at the end, \a inners contain the current state of the
                N-dimensional loop.
     */          
    template<class T_OBJECT, class T_PREDICATE>
      class NDimIterator
      {
        public:
          //! Type of inner iterators
          typedef T_OBJECT t_Object;
          //! Type of the predicates.
          typedef T_PREDICATE t_Predicate;

        protected:
          //! current value within each nested loop
          std::vector<t_Object> iterators;
          //! Start value of each nested loop
          std::vector<t_Object> start;
          //! End value of each nested loop
          std::vector<t_Object> end;
          //! Current iterator within the container of iterators (eg current loop)
          typename std::vector<t_Object> :: iterator current_iterator;
          //! Last iterator within the container of iterators
          typename std::vector<t_Object> :: iterator last_iterator;
          //! Predicate which ends the loop.
          t_Predicate predicate;


        public:
          //! Constructor
          NDimIterator(){};
          //! Destructor
          ~NDimIterator(){};

          //! Adds an inner loop to existing loops.
          void add( t_Object _begin, t_Object _end );
          //! Resets the iterator to start
          void reset();
          //! Clears. No loops. No nothing.
          void clear(){ start.clear(); end.clear(); iterators.clear(); }
          //! Erases loop number \a _i
          void erase( types::t_unsigned _i);
          //! Returns a constant reference to the component \a _i of the iterator.
          const t_Object& access( types::t_unsigned _i) const;
          //! Returns a reference to the component \a _i of the iterator.
          t_Object& access( types::t_unsigned _i);
          //! Increments the N-dimensional iterator.
          bool operator++() { return increment(); };
          //! Increments the N-dimensional iterator.
          bool increment();
          //! Returns true if the iterator is still valid.
          bool is_valid();
          //! Initializes the loop.
          void init_loop();
          //! Returns the current iterator for loop accessing internal iterators.
          t_Object& get_current_iterator() { return *current_iterator; }
          //! Increments the loop accessing internal iterators.
          bool next_iterator()  { return ( (++current_iterator) != last_iterator ); }
          //! Returns true if there are still internal iterators to access.
          bool loop_is_valid() { return ( current_iterator != last_iterator ); }

       protected:
#         ifdef LADA_DEBUG
           //! Debug.
           bool check_valid() const;
#         endif 

      };


      template<typename T_OBJECT, typename T_PREDICATE>
        inline void NDimIterator<T_OBJECT, T_PREDICATE> :: add( t_Object _begin, t_Object _end )
        { 
          start.push_back( _begin );
          end.push_back( _end );
          iterators.push_back( _begin );
        }
      template<typename T_OBJECT, typename T_PREDICATE>
        inline void NDimIterator<T_OBJECT, T_PREDICATE> :: reset()
        { 
          LADA_ASSERT( check_valid(), "Iterator is invalid.\n");
          std::copy( start.begin(), start.end(),
                     iterators.begin() );
        }
      template<typename T_OBJECT, typename T_PREDICATE>
        inline bool NDimIterator<T_OBJECT, T_PREDICATE> :: is_valid()
        {
          LADA_ASSERT( check_valid(), "Iterator is invalid.\n");
          return predicate( *iterators.begin(), *end.begin() ); 
        }
      template<typename T_OBJECT, typename T_PREDICATE>
        inline void NDimIterator<T_OBJECT, T_PREDICATE> :: init_loop()
        {
          LADA_ASSERT( check_valid(), "Iterator is invalid.\n");
          current_iterator = iterators.begin();
          last_iterator = iterators.end();
        }
      template<typename T_OBJECT, typename T_PREDICATE>
        const typename NDimIterator<T_OBJECT, T_PREDICATE> :: t_Object&
        NDimIterator<T_OBJECT, T_PREDICATE> :: access( types::t_unsigned _i ) const
        { 
          LADA_ASSERT( check_valid(), "Iterator is invalid.\n");
          LADA_ASSERT(_i >= iterators.size(), "Index out of range.\n");
        
          return ( iterators[_i] );
        }
      template<typename T_OBJECT, typename T_PREDICATE>
        typename NDimIterator<T_OBJECT, T_PREDICATE> :: t_Object&
          NDimIterator<T_OBJECT, T_PREDICATE> :: access( types::t_unsigned _i )
          { 
            LADA_ASSERT( check_valid(), "Iterator is invalid.\n");
            LADA_ASSERT(_i >= iterators.size(), "Index out of range.\n");
         
            return ( iterators[_i] );
          }
      template<typename T_OBJECT, typename T_PREDICATE>
      void NDimIterator<T_OBJECT, T_PREDICATE> :: erase( types::t_unsigned _i )
      { 
        LADA_ASSERT( check_valid(), "Iterator is invalid.\n");
        LADA_ASSERT(_i >= iterators.size(), "Index out of range.\n");

        start.erase( start.begin() + _i );
        end.erase( end.begin() + _i );
        iterators.erase( iterators.begin() + _i );
      }

      template<typename T_OBJECT, typename T_PREDICATE>
      bool NDimIterator<T_OBJECT, T_PREDICATE> :: increment()
      {
        LADA_ASSERT( check_valid(), "Iterator is invalid.\n");
        typename std::vector<t_Object> :: iterator ri_it;
        typename std::vector<t_Object> :: iterator ri_end;
        typename std::vector<t_Object> :: iterator ri_itstart;
        typename std::vector<t_Object> :: iterator ri_itend;

        ri_it = iterators.end() - 1;
        ri_end = iterators.begin();
        ri_itstart = start.end() - 1;
        ri_itend = end.end() - 1;

        while ( (ri_it - ri_end) >= 0 )
        {
          ++(*ri_it);
          if ( not predicate( *ri_it, *ri_itend ) )
          {
            if ( ri_it == ri_end )
              return false;
            *ri_it = *ri_itstart;
            --ri_it; --ri_itstart; --ri_itend;
          }
          else 
            break;
        }

        return true;
      }



#     ifdef LADA_DEBUG
        template<typename T_OBJECT, typename T_PREDICATE>
        bool NDimIterator<T_OBJECT, T_PREDICATE> :: check_valid() const
        {
          types::t_unsigned a = start.size();
          types::t_unsigned b = end.size();
          types::t_unsigned c = iterators.size();

          if ( a == 0 or b == 0 or c == 0 )
            return false;

          return ( a == b and b == c );
        }
#     endif


  } // namespace opt
} // namespace LaDa
#endif // _NDIM_ITERATOR_H_
