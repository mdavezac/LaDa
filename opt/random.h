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
    //! Flips a coin, retuns true or false.
    bool flip();
    //! Flips a coin, returns +1 or -1.
    types::t_real rflip();
    //! Stops the radom number generator
    void destroy();
    //! calls generator
    types::t_real rng();
  }
}


#endif
