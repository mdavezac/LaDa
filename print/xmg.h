//
//  Version: $Id$
//
#ifndef __PRINT_XMG_H_
#define __PRINT_XMG_H_

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

#include "operations.h"

namespace Print
{
  std::string make_commented_string( const std::string& );
  class Xmg
  {
    template< class T_TYPE > friend Xmg& operator<< ( Xmg &, const T_TYPE & );
    friend std::string make_commented_string( const std::string &_str );
    protected:
      enum t_operation { COMMENT, CLEAR, INDENT, UNINDENT, ADDTOLAST, REMOVELAST,
                         CLEARALL };
    public:
      const static t_operation comment;
      const static t_operation clear; 
      const static t_operation indent;
      const static t_operation unindent; 
      const static t_operation addtolast;
      const static t_operation removelast;
      const static t_operation clearall;
    protected:
      const static std::string comment_string;

    protected:
      types::t_unsigned indentation;
      std::string filename;
      std::ofstream file;
      std::list< std::string > line_list;
      std::ostringstream stream;
      bool is_empty;
      bool do_print;

    public:
      Xmg() : indentation(0), filename(""), is_empty(true), do_print(false) {};
      ~Xmg() { close(); }
      
      bool open();
      void close();

      void clear_all()
      {
        if ( not do_print )  return;
        line_list.clear(); 
      }
      bool is_open() 
      {
        return do_print and file.is_open(); 
      }
      void add_line( const std::string &_str )
      {
        if ( not do_print )  return;
        line_list.push_back(_str); 
      }
      void add_comment( const std::string &_str);
      void add_to_last( const std::string &_str);
      void remove_last()
      {
        if ( not do_print )  return;
        line_list.pop_back(); 
      }
      void init(const std::string &_f);
      std::string get_filename() const { return filename; }

    protected:
      void flushall();
      template<class T_TYPE>
      void to_current_stream ( const T_TYPE& _type );
      void do_indent()
        { for( types::t_unsigned i = 0; i < indentation; ++i) stream << "  "; }
      void add_to_last();
      void special_op( t_operation _op);


#ifdef _MPI
    public:
#ifndef _PRINT_ALL_PROCS
      void sync_filename() {}
      void sync_filename( std::string & ) {}
#else
      void sync_filename();
      void sync_filename( std::string &_filename );
#endif
#endif
  };

  extern Xmg xmg;

  template<class T_TYPE> inline void Xmg :: to_current_stream ( const T_TYPE &_type)
  {
      if ( not do_print )  return;
      stream << _type;
  }

  template<> inline void Xmg :: to_current_stream<Xmg::t_operation>(const Xmg::t_operation &_op)
  {
    special_op( _op );
  }
  template<> inline void Xmg :: to_current_stream<Print::t_Operation>(const Print::t_Operation &_op)
  {
    if ( _op == Print::ENDL  and stream.str().empty() ) return;
    switch ( _op )
    {
      case ENDL: line_list.push_back( stream.str() ); stream.str(""); break;
      case FLUSH: flushall(); break;
      default: Print::apply_ops( stream, _op ); break;
    }
  }
  template<> inline void Xmg::to_current_stream<Print::setw>(const Print::setw &_op) { _op(stream); }
  template<> inline void Xmg::to_current_stream<Print::setfill>(const Print::setfill &_op) { _op(stream); }
  template<> inline void Xmg::to_current_stream<Print::setprecision>(const Print::setprecision &_op)
    { _op(stream); }
   

  template<class T_TYPE> Xmg& operator<< ( Xmg& _this, const T_TYPE &_whatever)
  {
    if ( not _this.do_print )  return _this;
    _this.to_current_stream( _whatever );
    return _this;
  }

}

#endif // _PRINT_XMGRACE_H_
