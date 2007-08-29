//
//  Version: $Id$
//
#ifndef _PRINT_MANIP_H_
#define _PRINT_MANIP_H_

#include <string>

namespace Print
{

  std::string StripDir( std::string _string );
  std::string StripDir( const std::string &_dir, const std::string &_str );
  std::string StripEdges( std::string _string );
  std::string reformat_home( std::string _str );

}
#endif 
