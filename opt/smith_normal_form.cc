//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "smith_normal_form.h"

//! \cond
extern "C" void FC_FUNC(smithnormalform, SMITHNORMALFORM)
                       ( const int*, int*, int*, int* );
//! \endcond

namespace LaDa
{

  namespace opt
  {
    
    void smith_normal_form( atat::iMatrix3d& _S, atat::iMatrix3d & _L,
                            const atat::iMatrix3d& _M, atat::iMatrix3d &_R )
    {
      types::t_int s[9], l[9], m[9], r[9];
      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
        {
          m[ j * 3 + i ] = _M(i,j);
          s[ j * 3 + i ] = 0;
          l[ j * 3 + i ] = 0;
          r[ j * 3 + i ] = 0;
        }
      FC_FUNC(smithnormalform, SMITHNORMALFORM)( m, l, s, r );
      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
        {
          _S(i,j) = s[ j * 3 + i ];
          _L(i,j) = l[ j * 3 + i ];
          _R(i,j) = r[ j * 3 + i ];
        }
    }
   
  } // namespace Crystal

} // namespace LaDa
