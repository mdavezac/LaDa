#ifndef _DARWIN_PRINT_XMGRACE_H_
#define _DARWIN_PRINT_XMGRACE_H_

#include <opt/types.h>
#include <string>
#include <fstream>
#include <list>

namespace darwin
{
  class PrintXmg
  {
    protected:
      types::t_unsigned indentation;
      std::string filename;
      std::ofstream file;
      std::list< std::string > line_list;

    public:
      PrintXmg(const std::string &_f) : indentation(0), filename(_f) {};
      ~PrintXmg() { close(); }
      
      bool open();
      void close();
      void flush();
      void clear()
        { line_list.clear(); }
      bool is_open() 
        { return file.is_open(); }
      void add_line( const std::string &_str )
        { line_list.push_back(_str); }
      void add_comment( const std::string &_str);
      void add_to_last( const std::string &_str);
      void indent()
        { ++indentation; }
      void deindent()
        { if ( indentation ) --indentation; }
      void remove_last()
        { line_list.pop_back(); }
      void init(const std::string &_f);
  };

  extern PrintXmg printxmg;
}

#endif // _PRINT_XMGRACE_H_
