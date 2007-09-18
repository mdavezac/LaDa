//
//  Version: $Id$
//
#ifndef _MOLECULARITY_IMPL_H_
#define _MOLECULARITY_IMPL_H_

#include <algorithm>
#include <functional>
#include <ext/algorithm>
#ifndef __PGI
  #include<ext/functional>
  using __gnu_cxx::compose1;
#else
  #include<functional>
  using std::compose1;
#endif

#include "print/xmg.h"
#include "print/stdout.h"
#include "functors.h"

namespace Molecularity
{
  template<class T_R_IT, class T_K_IT>
  void fourrier_to_kspace( T_R_IT _rfirst, T_R_IT _rend,
                           T_K_IT _kfirst, T_K_IT _kend ) // sets kvector values from rspace values
  {
    const std::complex<types::t_real>
       imath(0, -2*3.1415926535897932384626433832795028841971693993751058208);
    
    for (; _kfirst != _kend; ++_kfirst)
    {
      _kfirst->type = std::complex<types::t_real>(0);
      for(T_R_IT i_r( _rfirst ); i_r != _rend; i_r += 2 )
      {
        _kfirst->type +=   exp( imath * ( i_r->pos[0] * _kfirst->pos[0] +
                                          i_r->pos[1] * _kfirst->pos[1] +
                                          i_r->pos[2] * _kfirst->pos[2] ) )
                         * i_r->type;
      }
    }
  }
  template<class T_R_IT, class T_K_IT, class T_O_IT >
  void fourrier_to_rspace( T_R_IT _rfirst, T_R_IT _rend,
                           T_K_IT _kfirst, T_K_IT _kend,
                           T_O_IT _rout ) // sets rvector values from kspace values
  {
    const std::complex<types::t_real>
       imath(0, 2*3.1415926535897932384626433832795028841971693993751058208);
    for (; _rfirst != _rend; _rfirst+=2, ++_rout)
    {
      *_rout = 0.0;
      for(T_K_IT i_k=_kfirst; i_k != _kend; ++i_k)
      {
        *_rout +=   exp( imath * ( _rfirst->pos[0] * i_k->pos[0] +
                                   _rfirst->pos[1] * i_k->pos[1] +
                                   _rfirst->pos[2] * i_k->pos[2] ) )
                  * i_k->type;
      }
    }
  }

} // namespace TwoSites
#endif // _TWOSITES_IMPL_H_
