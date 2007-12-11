//
//  Version: $Id$
//
#include <string>
#include <stdexcept>
#include <iostream>

#include <print/manip.h>

#include "physics.h"

namespace Physics 
{
  namespace Atomic
  {
    types::t_int ExtractSymbol( char *_char, std::string &_s )
    { 
      types::t_int n = 0;
      while( *_char == ' ' or *_char == '-' )  { ++_char; ++n;  }
      if ( *_char == '\n'  )
      {
        std::cerr << "Encountered end of line prior to Atomic Symbol\n";
        return -n;
      }
      if( Print::is_lowercase( *_char ) ) return -n;
      _s.clear();
      _s.push_back( *_char );
      ++n; 
      if(     Print::is_lowercase( *(++_char) )
          and (*_char != ' ' and *_char != '-' and *_char  != '\n' )  )
      {
        _s.push_back( *_char );
        ++n; ++_char;  
      }
      if ( Z( _s ) ) return n;

      std::ostringstream sstr;
      std::cerr << "Could not determine atomic symbol from " << _s << std::endl;
      return -n; 
    }

  } // namespace Atomic
}
