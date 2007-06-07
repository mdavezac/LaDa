#ifndef _PHYSICAL_CONSTANTS_
#define _PHYSICAL_CONSTANTS_

#include <string>
#include <opt/types.h>

namespace physics 
{
  types::t_real a0( const std::string &_str );

  namespace atoms
  {
    types::t_unsigned Z(const std::string &_str);
    std::string Symbol(types::t_unsigned _n );
    types::t_unsigned Charge(types::t_unsigned _n );
    types::t_unsigned Charge( const std::string &_str );
  } // namespace atoms
}


#endif
