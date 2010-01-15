//
//  Version: $Id$
//
#ifndef _OPT_RANDOM_H_
#define _OPT_RANDOM_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>

namespace LaDa
{
  namespace math
  {
    //! Starts the radom number generator
    types::t_unsigned create( types::t_unsigned _seed = 0 );
    //! Starts the radom number generator
    types::t_unsigned seed( types::t_unsigned _seed = 0 );
    //! Flips a coin, retuns true or false.
    bool flip();
    //! Flips a coin, returns +1 or -1.
    types::t_real rflip();
    //! calls generator
    types::t_real rng();
    //! Gets a random number from a range [\a _first, \a_last).
    types::t_unsigned range( types::t_unsigned _first, types::t_unsigned _last );
    //! Gets a random number from a range [\a _first, \a_last).
    inline types::t_unsigned range( types::t_unsigned _last )
      { return range(0, _last ); }
  }
} // namespace Lada

#endif
