#ifndef _OPT_ALGORITHMS_H_
#define _OPT_ALGORITHMS_H_

#include "LaDaConfig.h"

#include<misc/types.h>
#include<opt/debug.h>

namespace LaDa
{
  namespace opt
  {
    template< class T_IT1, class T_IT2, class __APP >
      void concurrent_loop( T_IT1 _first, T_IT1 _last, T_IT2 _second, __APP _app )
      {
        for(; _first != _last; ++_first, ++_second )
          _app( *_first, *_second );
      }
    template< class T_IT1, class __APP >
      void concurrent_loop( T_IT1 _first, T_IT1 _last, types::t_unsigned _second, __APP _app )
      {
        for(; _first != _last; ++_first, ++_second )
          _app( *_first, _second );
      }
    template< class T_IT1, class T_IT2, class T_IT3, class __APP >
      void concurrent_loop( T_IT1 _first, T_IT1 _last,
                            T_IT2 _second, T_IT3 _third, __APP _app )
      {
        for(; _first != _last; ++_first, ++_second, ++_third )
          _app( *_first, *_second, *_third );
      }

    //! From boost sandbox.
    template<typename I, typename O, typename Pred> 
      O copy_if ( I first, I last, O res, Pred p) 
      {
        while (first != last)
        {
          if (p(*first)) 
            *res++ = *first; 
          ++first; 
        } 
        return res; 
      } 
  } // end of opt namespace
} // namespace LaDa
#endif
