#ifndef _MY_TYPES_H_
#define _MY_TYPES_H_

#include "LaDaConfig.h"

#include<complex>
#include<limits>

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
    const t_real tolerance = 1.e-12; //!< all purpose tolerance global
    //! roundoff term for numerical noise crap.
    types::t_real const roundoff(5e3 * std::numeric_limits<types::t_real>::epsilon());
  }
} // namespace LaDa
#endif
