//
//  Version: $Id$
//
#ifndef _PHYSICAL_CONSTANTS_
#define _PHYSICAL_CONSTANTS_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>

#include <opt/types.h>

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


#endif
