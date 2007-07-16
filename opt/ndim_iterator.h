// Provides a multidimensional iterator for dynamic loops
// _Tp, the building block iterator, requires the incrementation ++
// _Predicate is a binary predicate taking _Tp on input 
//  and returning a bool
//
//  The last iterator on the stack is the fastest rolling iterator
//  in the loop
//
//  Usage:
//    Ndim_Iterator<t_Object, t_pred> iterator;
//    iterator.add(start, end)
//      where start(end) is the value of the iterator at start (end)
//    do {  } while( ++iterator )
//      this loop will break when the value of the first added
//      iterator is such that Predicate( first, values of first at
//      end) is false
//
//   iterators are pre-incremented with ++
//   each time an iterator and its end value returns false when
//   applied to _Predicate, it is reseted to its value at start,
//   and the previous iterator on the stack is incremented
//
// Functions for looping over the iterators are provided
// with init_loop, next_iterator, loop_is_valid, get_current_iterator
// Usage: for( iterator.init_loop(); iterator.loop_is_valid();
//             iterator.next_iterator() )
//          do something with iterator.get_currrent_iterator(); 
//
//
//
#ifndef _NDIM_ITERATOR_H_
#define _NDIM_ITERATOR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <vector>
#include <algorithm>

#include "opt/types.h"


#ifdef _DEBUG_ITERATORS_
  #define _DEBUG_ITERATOR_CHECK_ {if(!check_valid())\
         { std::cerr << "Invalid use of Ndim_Iterator" << std::endl;\
           exit(0) } } 
#else
  #define _DEBUG_ITERATOR_CHECK_ {}
#endif // _DEBUG_ITERATORS_


namespace opt
{

  template<class _Tp, class _Op>
    class Ndim_Iterator
    {
      std::vector<_Tp> iterators;
      std::vector<_Tp> start;
      std::vector<_Tp> end;
      typename std::vector<_Tp> :: iterator current_iterator;
      typename std::vector<_Tp> :: iterator last_iterator;
      _Op _predicate;


      public:
        Ndim_Iterator(){};
        ~Ndim_Iterator(){};

//       Ndim_Iterator( const _Op _pred ) : _predicate(_pred) {};

        void add( _Tp _begin, _Tp _end )
        { 
          start.push_back( _begin );
          end.push_back( _end );
          iterators.push_back( _begin );
        }
        void reset()
        { 
          _DEBUG_ITERATOR_CHECK_;
          std::copy( start.begin(), start.end(),
                     iterators.begin() );
        }
        void clear()
        {   
          start.clear(); end.clear(); iterators.clear();
        }
        void erase( types::t_unsigned i);
        const _Tp& access( types::t_unsigned i) const;
        _Tp& access( types::t_unsigned i);
        bool operator++() { return increment(); };
        bool increment();
        bool is_valid()
        {
          _DEBUG_ITERATOR_CHECK_;
          return _predicate( *iterators.begin(), *end.begin() ); 
        }
        void init_loop()
        {
          _DEBUG_ITERATOR_CHECK_;
          current_iterator = iterators.begin();
          last_iterator = iterators.end();
        }
        _Tp& get_current_iterator()
          { return *current_iterator; }
        bool next_iterator()
          { return ( (++current_iterator) != last_iterator ); }
        bool loop_is_valid()
          { return ( current_iterator != last_iterator ); }

     protected:
       #ifdef _DEBUG_ITERATORS_
         bool check_valid() const;
       #endif // _DEBUG_ITERATORS_

    };


    template<typename _Tp, typename _Op>
    const _Tp& Ndim_Iterator<_Tp, _Op> :: access( types::t_unsigned i ) const
    { 
      #ifdef _DEBUG_ITERATOR_
        _DEBUG_ITERATOR_CHECK_;
        if (i >= iterators.size() ) 
        {
          std::cerr << "Cannot access iterator: value beyond range "
                    << " _Tp& Ndim_Iterator :: access( types::t_unsigned i ) const"
                    << std::endl;
          exit(0);
        }
      #endif // _DEBUG_ITERATOR_

      return ( iterators[i] );
    }
    template<typename _Tp, typename _Op>
    _Tp& Ndim_Iterator<_Tp, _Op> :: access( types::t_unsigned i )
    { 
      #ifdef _DEBUG_ITERATOR_
        _DEBUG_ITERATOR_CHECK_;
        if (i >= iterators.size() ) 
        {
          std::cerr << "Cannot access iterator: value beyond range "
                    << " _Tp& Ndim_Iterator :: access( types::t_unsigned i ) const"
                    << std::endl;
          exit(0);
        }
      #endif // _DEBUG_ITERATOR_

      return ( iterators[i] );
    }
    template<typename _Tp, typename _Op>
    void Ndim_Iterator<_Tp, _Op> :: erase( types::t_unsigned i )
    { 
      #ifdef _DEBUG_ITERATOR_
        _DEBUG_ITERATOR_CHECK_;
        if (i >= iterators.size() ) 
        {
          std::cerr << "Cannot access iterator: value beyond range "
                    << " void Ndim_Iterator :: erase( types::t_unsigned i )"
                    << std::endl;
          exit(0);
        }
      #endif // _DEBUG_ITERATOR_

      start.erase( start.begin() + i );
      end.erase( end.begin() + i );
      iterators.erase( iterators.begin() + i );
    }

    template<typename _Tp, typename _Op>
    bool Ndim_Iterator<_Tp, _Op> :: increment()
    {
      _DEBUG_ITERATOR_CHECK_;
      typename std::vector<_Tp> :: iterator ri_it;
      typename std::vector<_Tp> :: iterator ri_end;
      typename std::vector<_Tp> :: iterator ri_itstart;
      typename std::vector<_Tp> :: iterator ri_itend;

      ri_it = iterators.end() - 1;
      ri_end = iterators.begin();
      ri_itstart = start.end() - 1;
      ri_itend = end.end() - 1;

      while ( (ri_it - ri_end) >= 0 )
      {
        ++(*ri_it);
        if ( not _predicate( *ri_it, *ri_itend ) )
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



    #ifdef _DEBUG_ITERATORS_
      template<typename _Tp, typename _Op>
      bool Ndim_Iterator<_Tp, _Op> :: check_valid() const
      {
        types::t_unsigned a = start.size();
        types::t_unsigned b = end.size();
        types::t_unsigned c = iterators.size();

        if ( a == 0 or b == 0 or c == 0 )
          return false;

        return ( a == b and b == c );
      }
    #endif // _DEBUG_ITERATORS_


} // namespace opt
#endif // _NDIM_ITERATOR_H_
