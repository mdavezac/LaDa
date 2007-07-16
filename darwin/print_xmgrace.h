#ifndef _DARWIN_PRINT_XMGRACE_H_
#define _DARWIN_PRINT_XMGRACE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <fstream>
#include <list>

#include "opt/types.h"
#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

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
      {
#ifdef _MPI
        if ( mpi::main.rank() != mpi::ROOT_NODE )
          return;
#endif 
        line_list.clear(); 
      }
      bool is_open() 
      {
#ifdef _MPI
        if ( mpi::main.rank() != mpi::ROOT_NODE )
          return false;
#endif 
        return file.is_open(); 
      }
      void add_line( const std::string &_str )
      {
#ifdef _MPI
        if ( mpi::main.rank() != mpi::ROOT_NODE )
          return;
#endif 
        line_list.push_back(_str); 
      }
      void add_comment( const std::string &_str);
      void add_to_last( const std::string &_str);
      void indent()
      {
#ifdef _MPI
        if ( mpi::main.rank() != mpi::ROOT_NODE )
          return;
#endif 
        ++indentation;
      }
      void deindent()
      {
#ifdef _MPI
        if ( mpi::main.rank() != mpi::ROOT_NODE )
          return;
#endif 
        if ( indentation ) --indentation; 
      }
      void remove_last()
      {
#ifdef _MPI
        if ( mpi::main.rank() != mpi::ROOT_NODE )
          return;
#endif 
        line_list.pop_back(); 
      }
      void init(const std::string &_f);
  };

  extern PrintXmg printxmg;
}

#endif // _PRINT_XMGRACE_H_
