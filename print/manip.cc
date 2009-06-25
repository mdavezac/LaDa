//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <cstdlib>
#include <sstream>

#include "manip.h"


namespace LaDa
{
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
      if ( _string.find_first_of(" \t\n") != std::string::npos )
        _string.erase(0, _string.find_first_not_of(" \t\n") );
      size_t l = _string.length();
      size_t t = _string.find_last_not_of(" \t\n") + 1;
      if ( t > l ) return _string;
      _string.erase(t, l );
      return _string;
    }

    std::string reformat_home ( const std::string &_str )
    {
      // Trim leading and tailing spaces
      std::string str = StripEdges(_str);
      if ( str[0] != '~' ) return str;
      std::string home="";
      if ( getenv("HOME") ) home = getenv("HOME");
      else if ( getenv("home") ) home = getenv("home");
      else return str;
      if ( str.find_first_of("/") )
        str.erase(0, str.find_first_of("/") );
      str.insert(0, home );
      return str;
    }

    std::string lowercase(std::string _string)
    {//change each element of the string to upper case
       for(unsigned int i=0;i<_string.length();i++)
          _string[i] = std::tolower(_string[i]);
       return _string;//return the converted string
    }
    std::string uppercase(std::string _string)
    {//change each element of the string to upper case
       for(unsigned int i=0;i<_string.length();i++)
          _string[i] = std::toupper(_string[i]);
       return _string;//return the converted string
    }

     
  // void BackwardFind( ifstream &_file, const std::string _comp )
  // {
  //   types::t_unsigned size = 10*_comp.size() * sizeof( char );
  //   if( size < 100 ) size = 100;
  //   types::t_int i;
  //
  //   // go to end of file (less eof char)
  //   _file.seekg(0, ios::end);
  //   i = (types::t_int) in.tellg(); // see how many bytes in file
  //   i -= size;                // backup before eof
  //   if ( i < 0 ) i = 0;
  //   while(
  // }

    std::string indent( const std::string &_ind, const std::string &_str )
    {
      if ( _str.empty() ) return _ind;
      std::ostringstream result;
      result << _ind;

      types::t_int old = 0;
      types::t_unsigned i = _str.find('\n', 0);
      while ( _str.find('\n', old ) != std::string::npos  )
      {
        result << _str.substr(old, i - old +1 ) << _ind;
        ++i; old = i;
        i = _str.find('\n', old);
      }
      result << _str.substr(old);
      return result.str();
    }

  }
} // namespace LaDa
