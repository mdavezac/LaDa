//
//  Version: $Id$
//
#include <stdlib.h>

#include "manip.h"

namespace Print
{
  std::string StripDir( std::string _string )
  {
    _string = StripEdges( _string );
    size_t t = _string.find_last_of("/");
    if ( t == std::string::npos )
      return _string;
    
    return _string.erase( 0, t + 1 );
  }
  std::string StripDir( const std::string &_dir, const std::string &_str )
  {
    std::string dir = reformat_home(_dir);
    std::string str = reformat_home(_str);
    if( str.find( dir ) == std::string::npos )
      return str;
    
    str.erase(0, dir.length() );
    if ( _str.find_first_of(" /\t\n") )
      str.erase(0, str.find_first_not_of(" /\t\n"));
    return str; 
  }
  std::string StripEdges( std::string _string )
  {
    if ( _string.find_first_of(" \t\n") )
      _string.erase(0, _string.find_first_not_of(" \t\n") );
    size_t l = _string.length();
    size_t t = _string.find_last_not_of(" \t\n") + 1;
    if ( t > l ) return _string;
    _string.erase(t, l );
    return _string;
  }

  std::string reformat_home ( std::string _str )
  {
    // Trim leading and tailing spaces
    _str = StripEdges(_str);
    if ( _str[0] != '~' ) return _str;
    std::string home="";
    if ( getenv("HOME") ) home = getenv("HOME");
    else if ( getenv("home") ) home = getenv("home");
    else return _str;
    if ( _str.find_first_of("/") )
      _str.erase(0, _str.find_first_of("/") );
    _str.insert(0, home );
    return _str;
  }



}

