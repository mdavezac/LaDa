//
//  Version: $Id$
//
#ifndef _PHYSICAL_CONSTANTS_
#define _PHYSICAL_CONSTANTS_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdexcept>

#include <string>

#include <opt/types.h>

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
  types::t_real planck( const std::string &_str );
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

    //! \brief Returns true if an atomic symbol is found in the next non-blank
    //!        and non '-' characters of \a _char 
    //! \details If '\n' is found before any non-blank or non-'-' character,
    //!          then the %function return false. If the second character found
    //!          is lowercase and is neither '-', nor ' ', nor '\n', then it is
    //!          assumed that the atomic symbol contains two characters. 
    //!          At this point \a _s contains the assumed atomic symbol. As a
    //!          last check, the function calls Physics::Atomic::Z(). If its
    //!          return is non-zero, ExtractSymbol() returns successfully.
    //! \return The funtions returns a negative number on failure, and a
    //!         positive number of success. This number indicates how far into
    //!         the string pointed to by \a _char the %function has gone. On
    //!         success, it should return one past the last letter of the
    //!         atomic symbol.
    types::t_int ExtractSymbol( char *_char, std::string &_s );
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

  inline types::t_real planck( const std::string &_str )
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
    inline std::string Symbol(types::t_unsigned _n )
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
    inline types::t_unsigned Charge(types::t_unsigned _n )
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
    inline types::t_unsigned Charge( const std::string &_str )
    {
      return Charge( Z( _str ) );
    }
  }
}
//!\endcond

#endif
