//
//  Version: $Id$
//
#ifndef _OPT_RANDOM_H_
#define _OPT_RANDOM_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "types.h"

namespace opt
{
  namespace random
  {
    //! Starts the radom number generator
    void create( types::t_unsigned _seed = 0 );
    //! Starts the radom number generator
    void seed( types::t_unsigned _seed = 0 );
    //! Stops the radom number generator
    void destroy();
    //! calls generator
    types::t_real rng();
  }
}


#endif
