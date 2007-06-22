#ifndef __MACHDEP_H__
#define __MACHDEP_H__

//#include <values.h>
#include <limits.h>
#include <float.h>
#define MAXINT INT_MAX
#define MAXFLOAT FLT_MAX

#include <opt/types.h>

namespace atat
{ 

#define Real double
#ifdef OLD_COMPLEX
  #define Complex complex
#else
  #define Complex complex<types::t_real>
#endif
#define PATHSEP '/'
#include <unistd.h>
#include <stdlib.h>

#include "fixagg.h"

using namespace std;


} // namespace atat

#endif
