//
//  Version: $Id$
//

#ifndef __PYTHONLADA_MISC_HPP_
#define __PYTHONLADA_MISC_HPP_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
#include <string>

namespace PythonLaDa
{
  template< class T_TYPE >
  std::string print( const T_TYPE &_at )
  { 
    std::ostringstream sstr;
    _at.print_out( sstr );
    return sstr.str();
  }
}

#endif // __PYTHONLADA_MISC_HPP_
