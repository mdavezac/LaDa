//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>

#include <eo/utils/eoRNG.h>

#include <opt/types.h>

#include "debug.h"

std::string keycheck( types::t_int _s)
{
  std::string result; result.resize(_s);
  for (types::t_int i=0 ; i < _s; ++i)
    result[i] = rng.flip() ? '1' : '0';
  return result;
}
