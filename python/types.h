#ifndef _MY_TYPES_H_
#define _MY_TYPES_H_

#include "LaDaConfig.h"

#include <complex>
#include <limits>

namespace LaDa
{
  //! Names a few simple variable types and globals for portability purposes
  namespace types 
  {
    typedef unsigned t_unsigned; //!< the unsigned integer type
    typedef int t_int;           //!< the signed integer type
    typedef double t_real;       //!< the real value type
    typedef char t_char;         //!< the character type, unused
    typedef std::complex<types::t_real> t_complex; //!< a complex real type
    //! \brief all-purpose global tolerance.
    //! \warning Setting this to a very small value may lead to bizarre errors.
    //!          For instance, some of the tests are known to fail occasionally
    //!          on a 64bit linux when tolerance = 1e-12.
    const t_real tolerance = 1.e-8; 
    //! roundoff term for numerical noise crap.
    types::t_real const roundoff(5e3 * std::numeric_limits<types::t_real>::epsilon());
  }
} // namespace LaDa
#endif
