#ifndef _PHYSICAL_CONSTANTS_
#define _PHYSICAL_CONSTANTS_

#include "LaDaConfig.h"

#ifdef __PGI
#define __BREAK 
#else
#define __BREAK break;
#endif

#include <stdexcept>

#include <string>

#include <opt/types.h>

namespace LaDa
{
  //! For mathematical constants
  namespace Math
  {
    //! \f$\pi\f$
    const types::t_real pi = 3.1415926535897932384626433832795028841971693993751058209749445920;
  }

  //! Physical constants and affiliate routines.
  namespace Physics 
  {
    //! Returns the the bhor radius \f$a_0\f$ in "A", "nm", "m", or "cm".
    types::t_real a0( const std::string &_str );
    //! Returns Hartree energy in "Ev", "Ry", and "H"
    types::t_real Hartree( const std::string &_str );
    //! Returns Rydberg energy in "Ev", "Ry", and "H"
    types::t_real Rydberg( const std::string &_str );
    //! Returns Planck's constant in "erg*s", "J*s", "eV*s", "Ry", and "H"
    types::t_real Planck( const std::string &_str );
    //! Returns \f$\hbar\f$ in "erg*s", "J*s", "eV*s", "Ry", and "H"
    types::t_real hbar( const std::string &_str );
    //! Returns electron volts in "eV", "erg", and "J"
    types::t_real eV( const std::string &_str );
    //! \brief Returns the mass of the electron \f$m_e\f$ in "kg", "g", "amu", "eV", and "MeV"
    //! \details "amu" stands for atomic mass unit. The conversion to eV is done
    //! via \f$m_e c^2 /\mathrm{"1 eV"}\f$
    types::t_real emass( const std::string &_str );

    
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
      //! Returns the atomic mass from the atomic symbol
      types::t_real Mass(const std::string &_str);

    } // namespace atoms
  }

  //! \cond
  namespace Physics 
  {
    inline types::t_real a0( const std::string &_str )
    {
      if ( _str == "A" )       return 0.529177249;
      else if ( _str == "nm" ) return 5.29177249;
      else if ( _str == "m" )  return 0.529177249e-8;
      else if ( _str == "cm" ) return 0.529177249e-6;
      else if ( _str == "au" or _str == "a.u." ) return 1e0;
      else throw std::runtime_error( "Unknown Unit for Bhor radius" );
    }

    inline types::t_real Hartree( const std::string &_str )
    {
      if ( _str == "eV" )      return 27.211396;
      else if ( _str == "Ry" ) return  2.0;
      else if ( _str == "H" )  return  1.0;
      else throw std::runtime_error( "Unknown Unit for Hartree constant" );
    }
    
    inline types::t_real Rydberg( const std::string &_str )
    {
      if ( _str == "eV" ) return 13.6056980;
      if ( _str == "Ry" ) return  1.0;
      if ( _str == "H" )  return  0.5;
      else throw std::runtime_error( "Unknown Unit for Rydberg constant" );
    }

    inline types::t_real Planck( const std::string &_str )
    {
      if( _str == "erg*s" ) return 6.626075510e-27;
      else if( _str == "J*s" ) return 6.626075510e-34;
      else if( _str == "eV*s" ) return 4.1356692e-15;
      else if( _str == "Ry" ) return 2.0 * Math::pi;
      else if( _str == "H" ) return 2.0 * Math::pi;
      else throw std::runtime_error( "Unknown Unit for Planck constant" );
    }

    inline types::t_real hbar( const std::string &_str )
    {
      if( _str == "erg*s" ) return 1.05457266e-27;
      else if( _str == "J*s" ) return 1.05457266e-34;
      else if( _str == "eV*s" ) return 6.5821220e-16;
      else if( _str == "Ry*s" ) return 1.0;
      else if( _str == "H*s" ) return 1.0;
      else throw std::runtime_error( "Unknown Unit for hbar constant" );
    }

    inline types::t_real eV( const std::string &_str )
    {
      if( _str == "eV" ) return 1;
      else if( _str == "J" ) return 1.60217733e-19;
      else if( _str == "erg" ) return 1.60217733e-12;
      else throw std::runtime_error( "Unknown Unit for eV constant" );
    }

    inline types::t_real emass( const std::string &_str )
    {
      if( _str == "kg" ) return 9.1093897e-31;
      else if( _str == "g" ) return 9.1093897e-34;
      else if( _str == "amu" ) return 5.48579903e-4;
      else if( _str == "eV" ) return 0.51099906e6;
      else if( _str == "MeV" ) return 0.51099906;
      else throw std::runtime_error( "Unknown Unit for emass constant" );
    }

    namespace Atomic
    {
      inline types::t_unsigned Z(const std::string &_str)
      {
        if( _str ==  "Li" ) return  3;
        else if( _str ==  "C" ) return  6; 
        else if( _str ==  "N" ) return  7; 
        else if( _str ==  "O" ) return  8; 
        else if( _str ==  "Na" ) return 11;
        else if( _str ==  "Mg" ) return 12;
        else if( _str ==  "Al" ) return 13;
        else if( _str ==  "Si" ) return 14;
        else if( _str ==  "P" ) return 15; 
        else if( _str ==  "S" ) return 16; 
        else if( _str ==  "Cl" ) return 17; 
        else if( _str ==  "K" ) return 19; 
        else if( _str ==  "Sc" ) return 21; 
        else if( _str ==  "Ni" ) return 28;
        else if( _str ==  "Cu" ) return 29;
        else if( _str ==  "Zn" ) return 30;
        else if( _str ==  "Ga" ) return 31;
        else if( _str ==  "Ge" ) return 32;
        else if( _str ==  "As" ) return 33;
        else if( _str ==  "Se" ) return 34;
        else if( _str ==  "Br" ) return 35;
        else if( _str ==  "Rb" ) return 37;
        else if( _str ==  "Cd" ) return 48;
        else if( _str ==  "In" ) return 49;
        else if( _str ==  "Sn" ) return 50;
        else if( _str ==  "Sb" ) return 51;
        else if( _str ==  "Te" ) return 52;
        else if( _str ==  "Cs" ) return 55;
        else if( _str ==  "Au" ) return 79;
        else if( _str ==  "Hg" ) return 80;
        return 0;
      }
      inline types::t_real Mass(const std::string &_str)
      {
        if     ( _str ==  "Li" ) return   6.941;
        else if( _str ==   "C" ) return  12.0107; 
        else if( _str ==   "N" ) return  14.00674; 
        else if( _str ==   "O" ) return  15.9994; 
        else if( _str ==  "Na" ) return  22.98976; 
        else if( _str ==  "Mg" ) return  24.305;
        else if( _str ==  "Al" ) return  26.9815;
        else if( _str ==  "Si" ) return  28.0855;
        else if( _str ==   "P" ) return  30.97376; 
        else if( _str ==   "S" ) return  32.066; 
        else if( _str ==  "Cl" ) return  35.453; 
        else if( _str ==   "K" ) return  39.0983;
        else if( _str ==  "Sc" ) return  44.956;
        else if( _str ==  "Ni" ) return  58.6934;
        else if( _str ==  "Cu" ) return  63.546;
        else if( _str ==  "Zn" ) return  65.39;
        else if( _str ==  "Ga" ) return  69.723;
        else if( _str ==  "Ge" ) return  72.61;
        else if( _str ==  "As" ) return  74.92160;
        else if( _str ==  "Se" ) return  78.96;
        else if( _str ==  "Br" ) return  79.904;
        else if( _str ==  "Rb" ) return  85.4678;
        else if( _str ==  "Cd" ) return 112.411;
        else if( _str ==  "In" ) return 114.818;
        else if( _str ==  "Sn" ) return 118.710;
        else if( _str ==  "Sb" ) return 121.760;
        else if( _str ==  "Te" ) return 127.60;
        else if( _str ==  "Cs" ) return 132.9054519;
        else if( _str ==  "Au" ) return 196.96655;
        else if( _str ==  "Hg" ) return 200.59;
        return 0;
      }
      inline std::string Symbol(types::t_unsigned _n )
      {
        switch( _n )
        {
          case   3: return "Li"; __BREAK
          case   6: return  "C"; __BREAK
          case   7: return  "N"; __BREAK
          case   8: return  "O"; __BREAK
          case  11: return "Na"; __BREAK
          case  12: return "Mg"; __BREAK
          case  13: return "Al"; __BREAK
          case  14: return "Si"; __BREAK
          case  15: return  "P"; __BREAK
          case  16: return  "S"; __BREAK
          case  17: return "Cl"; __BREAK
          case  19: return  "K"; __BREAK
          case  21: return "Sc"; __BREAK
          case  28: return "Ni"; __BREAK
          case  29: return "Cu"; __BREAK
          case  30: return "Zn"; __BREAK
          case  31: return "Ga"; __BREAK
          case  32: return "Ge"; __BREAK
          case  33: return "As"; __BREAK
          case  34: return "Se"; __BREAK
          case  35: return "Br"; __BREAK
          case  37: return "Rb"; __BREAK
          case  48: return "Cd"; __BREAK
          case  49: return "In"; __BREAK
          case  50: return "Sn"; __BREAK
          case  51: return "Sb"; __BREAK
          case  52: return "Te"; __BREAK
          case  55: return "Cs"; __BREAK
          case  79: return "Au"; __BREAK
          case  80: return "Hg"; __BREAK
        }
        return "error";
      }
      inline types::t_unsigned Charge(types::t_unsigned _n )
      {
        switch( _n )
        {
          case   3: return  1; __BREAK
          case   6: return  4; __BREAK
          case   7: return  5; __BREAK
          case   8: return  6; __BREAK
          case  12: return  2; __BREAK
          case  13: return  3; __BREAK
          case  14: return  4; __BREAK
          case  15: return  5; __BREAK
          case  16: return  6; __BREAK
          case  21: return  3; __BREAK
          case  28: return  9; __BREAK
          case  29: return 10; __BREAK
          case  30: return 11; __BREAK
          case  31: return  3; __BREAK
          case  32: return  4; __BREAK
          case  33: return  5; __BREAK
          case  34: return  6; __BREAK
          case  35: return  7; __BREAK
          case  48: return 11; __BREAK
          case  49: return  3; __BREAK
          case  50: return  4; __BREAK
          case  51: return  5; __BREAK
          case  52: return  6; __BREAK
          case  55: return  1; __BREAK
          case  79: return 10; __BREAK
          case  80: return 11; __BREAK
        }
        return 0;
      }
      inline types::t_unsigned Charge( const std::string &_str )
      {
        return Charge( Z( _str ) );
      }
    }
  }
  //!\endcond
} // namespace LaDa

#endif
