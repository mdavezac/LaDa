//
//  Version: $Id$
//
#ifndef _INTEGER__H_
#define _INTEGER__H_

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include <math.h>
#include "misc.h"

#include <opt/types.h>

namespace LaDa 
{ 
namespace atat
{ 

extern types::t_real zero_tolerance;

inline types::t_int near_zero(types::t_real x) {
  return (fabs(x)<zero_tolerance);
}

inline types::t_int is_int(types::t_real x) {
  return (fabs(x-rint(x)) < zero_tolerance);
}

types::t_int least_common_multiple(types::t_int a, types::t_int b);
types::t_int integer_ratio(types::t_int *p, types::t_int *q, types::t_real x, types::t_real epsilon);


} // namespace atat

} // namespace LaDa

#endif
