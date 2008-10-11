//
//  Version: $Id$
//
#ifndef __FIXAGG_H__
#define __FIXAGG_H__

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


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
