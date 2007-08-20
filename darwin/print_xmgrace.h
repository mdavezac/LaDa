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
    template< class T_TYPE > friend darwin::PrintXmg& operator<< ( darwin::PrintXmg &, const T_TYPE & );
    protected:
      enum t_operation { ENDL, COMMENT, CLEAR, FLUSH, INDENT, UNINDENT, ADDTOLAST, REMOVELAST,
                         CLEARALL };
    public:
      const static t_operation endl;
      const static t_operation comment;
      const static t_operation clear; 
      const static t_operation flush; 
      const static t_operation indent;
      const static t_operation unindent; 
      const static t_operation addtolast;
      const static t_operation removelast;
      const static t_operation clearall;
    protected:
      types::t_unsigned indentation;
      std::string filename;
      std::ofstream file;
      std::list< std::string > line_list;
      std::ostringstream stream;

    public:
      PrintXmg(const std::string &_f) : indentation(0), filename(_f) {};
      ~PrintXmg() { close(); }
      
      bool open();
      void close();

      void clear_all()
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
      void remove_last()
      {
#ifdef _MPI
        if ( mpi::main.rank() != mpi::ROOT_NODE )
          return;
#endif 
        line_list.pop_back(); 
      }
      void init(const std::string &_f);

    protected:
      void flushall();
      template<class T_TYPE>
      void to_current_stream ( const T_TYPE& _type );
      void do_indent()
        { for( types::t_unsigned i = 0; i < indentation; ++i) stream << "  "; }
      void add_to_last();
      void special_op( t_operation _op);
  };

  extern PrintXmg printxmg;

  template<class T_TYPE> inline void PrintXmg :: to_current_stream ( const T_TYPE &_type)
  {
#ifdef _MPI
    if ( not mpi::main.is_root_node() ) return;
#endif
      stream << _type;
  }

  template<> inline void PrintXmg :: to_current_stream<PrintXmg::t_operation>(const PrintXmg::t_operation &_op)
  {
    if ( _op == ENDL  and stream.str().empty() ) return;
    special_op( _op );
  }


  template<class T_TYPE> darwin::PrintXmg& operator<< ( darwin::PrintXmg& _this, const T_TYPE &_whatever)
  {
#ifdef _MPI
    if ( not mpi::main.is_root_node() ) return _this;
#endif
      _this.to_current_stream( _whatever );
    return _this;
  }
}


#endif // _PRINT_XMGRACE_H_
