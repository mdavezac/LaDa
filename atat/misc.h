//
//  Version: $Id$
//
#ifndef __MISC_H__
#define __MISC_H__

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <iostream>
#include <memory.h>
#include <sys/time.h>
#include "machdep.h"

#include <opt/types.h>

namespace atat
{ 

// -------------------------- divers

#define MAX_LINE_LEN 1024

#define SMOD(x,m) (((x)+(m)) % (m))
					// modulo signe, ex: SMOD(-1,4)=3

#define MAX(x,y) ((x)>(y) ? (x) : (y))
#define MIN(x,y) ((x)<(y) ? (x) : (y))
#define ABS(x) ((x)<0 ? (-(x)) : (x) )

template<class T>
inline T sqr(const T x) {return x*x;}

/*
template<class T>
T min(T a, T b) {
  return (a<b ? a : b);
}

template<class T>
T max(T a, T b) {
  return (a>b ? a : b);
}
*/

template<class T>
T sgn(T b) {
  return (b>(T)0 ? (T)1 : (b<(T)0 ? (T)(-1) : (T)0));
}

template<class T>
inline T ipow(T base, types::t_int exponent) {
  T accu=(T)1;
  while (exponent>0) {
    accu*=base;
    exponent--;
  }
  return accu;
}

template<class T>
inline void swap(T *a,T *b) {
  T tmp=*a;
  *a=*b;
  *b=tmp;
}

template<class T>
types::t_int is_between(T bw, T a, T b) {
  return (a<b ? (a<=bw && bw<=b) : (b<=bw && bw<=a));
}

#define countof(array) (sizeof(array)/sizeof(*array))

typedef short int shint;
typedef long int lint;
typedef unsigned int uint;
typedef unsigned long int ulint;
typedef unsigned char uchar;

#define MEMCPY(dest,src,cnt) memcpy(dest,src,(cnt)*sizeof(dest[0]))
#define MEMCLR(dest,cnt) memset(dest,0,(cnt)*sizeof(dest[0]))
#define MEMSET(dest,val,cnt) memset(dest,val,(cnt)*sizeof(dest[0]))
#define READARRAY(buf,len) read((char *)buf,len*sizeof(buf[0]))
#define WRITEARRAY(buf,len) write((char *)buf,len*sizeof(buf[0]))
#define READBINARY(var) read((char *)(&var),sizeof(var))
#define WRITEBINARY(var) write((char *)(&var),sizeof(var))

#define SWAPPTR(pt1,pt2) {void *tmp=(void *)(pt1); (pt1)=(pt2); (void *)(pt2)=tmp;}
#define ERRORQUIT(s) {cerr << s << endl; exit(1);}
#define COREDUMP {types::t_int *ptr=NULL; *ptr=0;}
#define ERRORQUITDUMP(s) {cerr << s << endl; COREDUMP;}

#define THIS (*this)

inline void rndseed(types::t_int seed=0) {
  if (seed==0) {
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv,&tz);
    srand(tv.tv_usec);
  }
  else {
    srand(seed);
  }
}

inline types::t_int random(types::t_int max) {
  return rand()%max;
}

inline types::t_real uniform01(void) {
  return (types::t_real)(rand())/(types::t_real)RAND_MAX;
}


} // namespace atat

#endif
