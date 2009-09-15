//
//  Version: $Id$
//
#ifndef LADA_ENUM_EXCEPTIONS_H_
#define LADA_ENUM_EXCEPTIONS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include <boost/throw_exception.hpp>
#include <boost/exception/error_info.hpp>
#include <boost/exception/info.hpp>


namespace LaDa
{
  namespace enumeration
  {
    //! Throws this when the lattice is not invariant through a given symmetry.
    struct internal: virtual boost::exception, virtual std::exception {}; 

    //! Throws this when the lattice is not invariant through a given symmetry.
    struct symmetry_not_of_lattice: virtual boost::exception, virtual std::exception {}; 
    //! Qualifies symmetry_not_of_lattice exceptions.
    typedef boost::error_info<struct enumeration,std::string> error_string; //(1)

    //! \brief Throws this when a transformation does not find equivalent supercell site.
    //! \details Most likely an internal error.
    struct symmetry_not_of_supercell : virtual symmetry_not_of_lattice {};
    
    //! \brief Thrown when the supercell under consideration is too large to be
    //!        represented by the integers in use.
    struct supercell_too_large : virtual boost::exception, virtual std::exception {}; 

    //! \brief Thrown when an integer is too large.
    struct integer_too_large : virtual boost::exception, virtual std::exception {}; 

    //! Thrown when the database is of incorrect size.
    struct incorrect_database : virtual boost::exception, virtual std::exception {}; 
    
  }
}

#endif
