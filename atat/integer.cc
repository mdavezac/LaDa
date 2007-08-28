//
//  Version: $Id$
//
#include "integer.h"


#include <opt/types.h>

namespace atat
{ 
types::t_real zero_tolerance=1e-3;

types::t_int least_common_multiple(types::t_int a, types::t_int b) {
  types::t_int test=2;
  types::t_int divider=1;
  while (test<=a && test<=b) {
    if ((a/test)*test == a && (b/test)*test == b) divider=test;
    test=test+1;
  }
  return a*b/divider;
}

types::t_real sign(types::t_real x) {
  return (x<0 ? -1. : (x>0 ? 1. : 0.));
}

types::t_int integer_ratio(types::t_int *pp, types::t_int *pq, types::t_real x, types::t_real epsilon) {
  types::t_int &p=*pp;
  types::t_int &q=*pq;
  types::t_int maxit=20;
  types::t_int p1,q1,p2,q2,n,g;
  types::t_real r;

  n=(types::t_int)rint(x);
  r=x-n;
  g=(types::t_int)sign(r);
  p=n;
  p1=p;
  q1=1;
  q=1;
  p2=1;
  q2=0;

  types::t_int counter=0;
  while (fabs((types::t_real)p/(types::t_real)q-x) > epsilon) {
    n=(types::t_int)rint(g/r);
    r=g/r-n;
    p=n*p1+g*p2;
    q=n*q1+g*q2;
    g=
    p2=p1;
    q2=q1;
    p1=p;
    q1=q;
    counter++;
    if (counter>maxit) return 0;
  }
  return 1;


} // namespace atat
}
