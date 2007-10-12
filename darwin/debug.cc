//
//  Version: $Id$
//
#ifndef _DEBUG_H_
#define _DEBUG_H_

#include <string>

#include <eo/utils/eoRNG.h>

#include <opt/types.h>

#include "debug.h"

std::string keycheck( types::t_int _s=20 )
{
  std::string result; result.resize(_s);
  for (types::t_int i=0 ; i < _s; ++i)
    result[i] = rng.flip() ? '1' : '0';
  return result;
}
#endif
