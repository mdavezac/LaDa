#ifndef _MY_TYPES_H_
#define _MY_TYPES_H_

#include<complex>

namespace types 
{
  typedef unsigned t_unsigned;
  typedef int t_int;
  typedef double t_real;
  typedef char t_char;
  const t_real tolerance = 0.000001;
  typedef std::complex<types::t_real> t_complex;
}
#endif
