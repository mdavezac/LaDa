//
//  Version: $Id$
//
#ifndef _LAPACK_C_H_
#define _LAPACK_C_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>
#include <atat/vectmac.h>

//! Lapack wrappers for C
namespace Lapack
{
  //! \brief Eigenvalues and Eigenvectors of symmetric 3d matrices
  //! \param[in] _mat matrix for which to find eigenvalues
  //! \param[out] _vects Eigenvectors
  //! \param[out] _vals Eigenvalues
  bool eigen( const atat::rMatrix3d &_mat, atat::rMatrix3d &_vecs, types::t_real _vals[3] );

}

#endif
