//
//  Version: $Id$
//
#ifndef _PHYSICAL_CONSTANTS_
#define _PHYSICAL_CONSTANTS_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>

#include "opt/types.h"

//! Physical constants and affiliate routines.
namespace Physics 
{
  //! Returns the the bhor radius in "A", "nm", "m", or "cm".
  types::t_real a0( const std::string &_str );
  
  //! Atomic numbers, symbols, and affiliate routines.
  namespace Atomic
  {
    //! Returns the atomic number from an atomic symbol
    types::t_unsigned Z(const std::string &_str);
    //! Returns an atomic symbol from an atomic number
    std::string Symbol(types::t_unsigned _n );
    //! Returns the number of valence electrons from the atomic number
    types::t_unsigned Charge(types::t_unsigned _n );
    //! Returns the number of valence electrons from the atomic symbol
    types::t_unsigned Charge( const std::string &_str );
  } // namespace atoms
}


#endif
