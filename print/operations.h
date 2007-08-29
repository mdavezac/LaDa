//
//  Version: $Id$
//
#ifndef __PRINT_OPERATIONS_H_
#define __PRINT_OPERATIONS_H_

#include <iomanip> 

namespace Print
{
  enum t_Operation { ENDL, FLUSH, LEFT, RIGHT, FIXED, SCIENTIFIC, BOOLALPHA, NOBOOLALPHA };
  const t_Operation endl        = ENDL;
  const t_Operation flush       = FLUSH;
  const t_Operation left        = LEFT;
  const t_Operation right       = RIGHT;
  const t_Operation fixed       = FIXED;
  const t_Operation scientific  = SCIENTIFIC;
  const t_Operation boolalpha   = BOOLALPHA;
  const t_Operation noboolalpha = NOBOOLALPHA;

  template <class T_TYPE> inline void apply_ops( T_TYPE &_type, const t_Operation &_op)
  {
    switch( _op )
    {
      case ENDL:        _type << std::endl;        break;
      case FLUSH:       _type << std::flush;       break;
      case LEFT:        _type << std::left;        break;
      case RIGHT:       _type << std::right;       break;
      case FIXED:       _type << std::fixed;       break;
      case SCIENTIFIC:  _type << std::scientific;  break;
      case BOOLALPHA:   _type << std::boolalpha;   break;
      case NOBOOLALPHA: _type << std::noboolalpha; break;
    }
  }

  class setw
  {
    public:
      types::t_unsigned width;
    public:
      setw(types::t_unsigned _w) : width(_w) {} 
      template<class T_TYPE> void operator()( T_TYPE &_whatever ) const
        { _whatever << std::setw(width); }
  };
  typedef setw width;
  class setprecision
  {
    public:
      types::t_unsigned prec;
    public:
      setprecision(types::t_unsigned _n) : prec(_n) {} 
      template<class T_TYPE> void operator()( T_TYPE &_whatever ) const
        { _whatever << std::setprecision(prec); }
  };
  class setfill
  {
    public:
      char c;
    public:
      setfill(char _c) : c(_c) {} 
      template<class T_TYPE> void operator()( T_TYPE &_whatever ) const
        { _whatever << std::setfill(c); }
  };
}

#endif
