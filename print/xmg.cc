//
//  Version: $Id$
//
#include<iostream>
#include <sstream>

#include <revision.h>

#include "xmg.h"
#include "manip.h"

namespace Print
{
  const Xmg :: t_operation Xmg :: comment    = Xmg :: COMMENT;
  const Xmg :: t_operation Xmg :: clear      = Xmg :: CLEAR;
  const Xmg :: t_operation Xmg :: indent     = Xmg :: INDENT;
  const Xmg :: t_operation Xmg :: unindent   = Xmg :: UNINDENT;
  const Xmg :: t_operation Xmg :: addtolast  = Xmg :: ADDTOLAST;
  const Xmg :: t_operation Xmg :: removelast = Xmg :: REMOVELAST;
  const Xmg :: t_operation Xmg :: clearall   = Xmg :: CLEARALL;
  const std::string Xmg :: comment_string = "# ";

  void Xmg :: init (const std::string &_f)
  { 
#ifdef _MPI
    if ( not mpi::main.is_root_node() ) return;
#endif 
    filename = reformat_home( _f );
    file.open( filename.c_str(), std::ios_base::out|std::ios_base::trunc ); 
    is_empty = true;
    if (file.fail() )
    {
      std::cerr << "Could not open " << filename << std::endl;
      return;
    }
    do_print = true;
    close();
  }
  bool Xmg :: open ()
  {
    if ( not do_print ) return true;
    if ( file.is_open() ) return true;
    file.open( filename.c_str(), std::ios_base::out|std::ios_base::app ); 
    return file.is_open();
  }    
  void Xmg :: close ()
  {
    if ( not do_print ) return;
    if ( not file.is_open() ) return;
    flushall();
    file.flush();
    file.close();
  }    
  void Xmg :: flushall ()
  {
    if ( not do_print ) return;
    if ( line_list.empty() ) return; 

    open();

    if ( is_empty )  // print revision
      file << comment_string << "Subversion Revision " << SVN::Revision << std::endl;

    is_empty = false;

    std::list< std::string > :: const_iterator i_str = line_list.begin();
    std::list< std::string > :: const_iterator i_str_end = line_list.end();
    for(; i_str != i_str_end; ++i_str) 
      file << *i_str << "\n";

    line_list.clear();
    file.flush();
  }    

  void Xmg :: add_comment( const std::string &_str )
  {
    if ( not do_print ) return;
    std::ostringstream sstr;
    sstr << comment_string;
    for( types::t_unsigned i = 0; i < indentation; ++i)
      sstr << "  ";
    sstr << _str;
    std::string str = sstr.str();
    line_list.push_back( str );
  }
  void Xmg :: add_to_last( const std::string &_str )
  {
    if ( not do_print ) return;
    if( line_list.empty() ) return;
    std::string &str = line_list.back();
    str += _str;
  }
  void Xmg :: add_to_last()
  {
    if ( not do_print ) return;
    if( line_list.empty() ) return;

    stream.str(""); stream << line_list.back();
    line_list.pop_back();
  }

  void Xmg :: special_op( t_operation _op )
  {
    switch ( _op )
    {
      case CLEAR: stream.str(""); break;
      case COMMENT: stream.str(""); stream << comment_string; do_indent(); break;
      case INDENT: ++indentation; break;
      case UNINDENT: --indentation; break;
      case ADDTOLAST: add_to_last(); break;
      case REMOVELAST: line_list.pop_back(); break;
      case CLEARALL: stream.str(""); clear_all(); break;
    }
  }


#ifdef _MPI
#ifdef _PRINT_ALL_PROCS
  void Xmg::sync_filename( std::string &_filename )
  {
    filename = _filename;
    sync_filename();
  }
  void Xmg::sync_filename()
  {
    mpi::BroadCast bc( mpi::main );
    bc << filename
       << mpi::BroadCast::allocate
       << filename
       << mpi::BroadCast::broadcast
       << filename;

    if ( mpi::main.is_root_node() ) return;

    init(filename);
  }
#endif
#endif



  std::string make_commented_string( const std::string &_str )
  {
    std::string copy(_str);
    std::ostringstream result;
    if ( copy.empty() ) return "";
    result << Xmg::comment_string;

    types::t_int old = 0;
    types::t_int i = copy.find('\n', 0);
    while ( i != std::string::npos  )
    {
      result << copy.substr(old, i - old +1 ) << Xmg::comment_string;
      ++i; old = i;
      i = copy.find('\n', old);
    }
    result << copy.substr(old);
    return result.str();
  }

  Xmg xmg;
}
