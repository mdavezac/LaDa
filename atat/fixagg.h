//
//  Version: $Id$
//
#ifndef __FIXAGG_H__
#define __FIXAGG_H__

#include <opt/types.h>

namespace atat
{ 

template<class T>
class Aggregate {
public:
  T x;
  Aggregate(void) {}
  Aggregate(const T &r) {x=r;}
  operator T & () {return x;}
};


} // namespace atat

#endif
