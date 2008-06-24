//
//  Version: $Id$
//
#ifndef __VECTMAC_H__
#define __VECTMAC_H__

#define STREAM_VECTOR
#include "fxvector.h"

#include <opt/types.h>

namespace atat
{ 

typedef Vector2d<types::t_real> rVector2d;
typedef Vector2d<types::t_int> iVector2d;
typedef Vector2d<types::t_real> rVector2d;
typedef Vector2d<types::t_int> iVector2d;
typedef Matrix2d<types::t_real> rMatrix2d;
typedef Matrix2d<types::t_int> iMatrix2d;
typedef BoundingBox<types::t_int,2> iBoundingBox2d;

typedef Vector3d<types::t_real> rVector3d;
typedef Vector3d<types::t_int> iVector3d;
typedef Vector3d<types::t_real> rVector3d;
typedef Vector3d<types::t_int> iVector3d;
typedef Matrix3d<types::t_real> rMatrix3d;
typedef Matrix3d<types::t_int> iMatrix3d;
typedef BoundingBox<types::t_int,3> iBoundingBox3d;
typedef BoundingBox<types::t_real,3> rBoundingBox3d;

inline rVector2d to_real(const iVector2d &v) {
  return rVector2d((types::t_real)v(0),(types::t_real)v(1));
}

inline rVector3d to_real(const iVector3d &v) {
  return rVector3d((types::t_real)v(0),(types::t_real)v(1),(types::t_real)v(2));
}

inline rMatrix3d to_real(const iMatrix3d &a) {
  rMatrix3d b;
  for (types::t_int i=0; i<3; i++) {
    for (types::t_int j=0; j<3; j++) {
      b(i,j)=(types::t_real)a(i,j);
    }
  }
  return rMatrix3d(b);
}

inline types::t_int iround(types::t_real x) {
  return (types::t_int)rint(x);
}

inline iVector2d to_int(const rVector2d &v) {
  return iVector2d(iround(v(0)),iround(v(1)));
}

inline iVector3d to_int(const rVector3d &v) {
  return iVector3d(iround(v(0)),iround(v(1)),iround(v(2)));
}

inline iMatrix3d to_int(const rMatrix3d &a) {
  iMatrix3d b;
  for (types::t_int i=0; i<3; i++) {
    for (types::t_int j=0; j<3; j++) {
      b(i,j)=iround(a(i,j));
    }
  }
  return b;
}


} // namespace atat

// includes return type declaration for usual operators.
#include "lambda.impl.h"

#endif
