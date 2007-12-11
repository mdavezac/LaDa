//
//  Version: $Id$
//

#include "lapack.h"

#include <opt/traits.h>

//! \cond
extern "C" void FC_FUNC(dsyev, DSYEV)( char*, char*, int*, types::t_real*, 
                                       int*, types::t_real*, types::t_real*, 
                                       int*, int* );
//! \endcond

namespace Lapack
{
  bool eigen( const atat::rMatrix3d &_mat, atat::rMatrix3d &_vecs, types::t_real _vals[3] )
  {
    double work[9];
    double inv[9];
    double eig[3];
    
    int info;
    char job('V'), mattype('U');
    int three(3), nine(9);
    _vecs = _mat;

    for( types::t_unsigned i=0; i<3; ++i )
      for( types::t_unsigned j=0; j<3; ++j )
        if(  opt::Fuzzy<types::t_real>::equal( _mat(i,j), 0.0 ) ) inv[3*i+j] = 0.0;
        else inv[3*i+j] = (double) _mat(i,j);

    FC_FUNC(dsyev, DSYEV)( &job, &mattype, &three, inv, &three, eig, work, &nine, &info);

    for( types::t_unsigned i=0; i<3; ++i )
    {
      for( types::t_unsigned j=0; j<3; ++j )
        if( opt::Fuzzy<double>::equal( inv[3*i+j], 0.0 ) ) _vecs(i,j) = 0.0;
        else _vecs(i,j) = (types::t_real) inv[3*i+j];

      if( opt::Fuzzy<double>::equal( eig[i], 0.0 ) ) _vals[i] = 0.0;
      else _vals[i] = (types::t_real) eig[i];
    }
    

    return info == 0;
  }
}
