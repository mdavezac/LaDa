#ifndef _OPT_MINIMIZE_H_
#define _OPT_MINIMIZE_H_

// includes all minimizers
#include "opt/opt_minimize_linear.h"
// #include <opt/opt_minimize_wang.h>
// #include <opt/opt_minimize_ssquared.h>
#include "opt/types.h"

namespace opt
{
  const types::t_unsigned NO_MINIMIZER       = 0;
  const types::t_unsigned WANG_MINIMIZER     = 1;
  const types::t_unsigned PHYSICAL_MINIMIZER = 2;
  const types::t_unsigned LINEAR_MINIMIZER   = 3;
  const types::t_unsigned SA_MINIMIZER       = 4;
}

#endif
