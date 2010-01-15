//
//  Version: $Id$
//
#ifndef LADA_MATH_IPOW_H
#define LADA_MATH_IPOW_H

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <opt/types.h>

namespace LaDa
{

  namespace math
  { 

    //! Exponentiation for integer powers.
    template<class T>
      inline T ipow(T base, size_t exponent) 
      {
        T result(1);
        for(; exponent; --exponent, result *= base);
        return result;
      }

  } // namespace math

} // namespace LaDa

#endif
