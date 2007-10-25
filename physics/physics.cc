//
//  Version: $Id$
//
#ifndef _PHYSICAL_CONSTANTS_
#define _PHYSICAL_CONSTANTS_

#include <string>

#include "opt/types.h"
#include "physics.h"

namespace Physics 
{
  types::t_real a0( const std::string &_str )
  {
    if ( _str == "A" )
      return 0.529177249;
    if ( _str == "nm" )
      return 5.29177249;
    else if ( _str == "m" )
      return 0.529177249e-8;
    else if ( _str == "cm" )
      return 0.529177249e-6;
    return 1.0;
  }
  
  namespace Atomic
  {
    types::t_unsigned Z(const std::string &_str)
    {
      if( _str ==  "Li" ) return  3;
      else if( _str ==  "C" ) return  6; 
      else if( _str ==  "N" ) return  7; 
      else if( _str ==  "O" ) return  8; 
      else if( _str ==  "Al" ) return 13;
      else if( _str ==  "Si" ) return 14;
      else if( _str ==  "P" ) return 15; 
      else if( _str ==  "S" ) return 16; 
      else if( _str ==  "Ni" ) return 28;
      else if( _str ==  "Cu" ) return 29;
      else if( _str ==  "Zn" ) return 30;
      else if( _str ==  "Ga" ) return 31;
      else if( _str ==  "Ge" ) return 32;
      else if( _str ==  "As" ) return 33;
      else if( _str ==  "Se" ) return 34;
      else if( _str ==  "Br" ) return 35;
      else if( _str ==  "Cd" ) return 48;
      else if( _str ==  "In" ) return 49;
      else if( _str ==  "Sn" ) return 50;
      else if( _str ==  "Sb" ) return 51;
      else if( _str ==  "Te" ) return 52;
      else if( _str ==  "Au" ) return 79;
      else if( _str ==  "Hg" ) return 80;
      return 0;
    }
    std::string Symbol(types::t_unsigned _n )
    {
      switch( _n )
      {
        case   3: return "Li"; break;
        case   6:  return "C"; break;
        case   7:  return "N"; break;
        case   8:  return "O"; break;
        case  13: return "Al"; break;
        case  14: return "Si"; break;
        case  15:  return "P"; break;
        case  16:  return "S"; break;
        case  28: return "Ni"; break;
        case  29: return "Cu"; break;
        case  30: return "Zn"; break;
        case  31: return "Ga"; break;
        case  32: return "Ge"; break;
        case  33: return "As"; break;
        case  34: return "Se"; break;
        case  35: return "Br"; break;
        case  48: return "Cd"; break;
        case  49: return "In"; break;
        case  50: return "Sn"; break;
        case  51: return "Sb"; break;
        case  52: return "Te"; break;
        case  79: return "Au"; break;
        case  80: return "Hg"; break;
      }
      return "error";
    }
    types::t_unsigned Charge(types::t_unsigned _n )
    {
      switch( _n )
      {
        case   3: return  1; break;
        case   6: return  4; break;
        case   7: return  5; break;
        case   8: return  6; break;
        case  13: return  3; break;
        case  14: return  4; break;
        case  15: return  5; break;
        case  16: return  6; break;
        case  28: return  9; break;
        case  29: return 10; break;
        case  30: return 11; break;
        case  31: return  3; break;
        case  32: return  4; break;
        case  33: return  5; break;
        case  34: return  6; break;
        case  35: return  7; break;
        case  48: return 11; break;
        case  49: return  3; break;
        case  50: return  4; break;
        case  51: return  5; break;
        case  52: return  6; break;
        case  79: return 10; break;
        case  80: return 11; break;
      }
      return 0;
    }
    types::t_unsigned Charge( const std::string &_str )
    {
      return Charge( Z( _str ) );
    }
  } // namespace Atomic
}


#endif
