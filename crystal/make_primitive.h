//
//  Version: $Id$
//
#ifndef LADA_CRYSTAL_MAKE_PRIMITIVE_H_
#define LADA_CRYSTAL_MAKE_PRIMITIVE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/shared_ptr.hpp>

#include "lattice.h"

namespace LaDa
{
  namespace Crystal 
  {
    //! Creates a primitive lattice.
    boost::shared_ptr<Lattice> make_primitive( Lattice const &_lattice );

  } // namespace Crystal
} // namespace LaDa
#endif
