#include "LaDaConfig.h"

#ifndef _CUBIC_CE_
# error Undefined _CUBIC_CE_
#endif

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/errors.hpp>

#include <opt/types.h>
#include <opt/debug.h>

#include "ce.hpp"
#include "ce.impl.hpp"

namespace LaDa
{
  namespace Python
  {

    void expose_ce()
    {
      detailsCE::expose_ce_functional< CE::ConstituentStrain::Harmonic::Cubic >
        ( "Cubic", "CE for cubic lattices.");
//     detailsCE::expose_ce_functional< CE::ConstituentStrain::Harmonic::Tetragonal >
//       ( "Tetra", "CE for tetragonal lattices.");
    }
  } // end of PythonLaDa namespace
} // namespace LaDa
