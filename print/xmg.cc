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
  const Xmg :: t_operation Xmg :: endl       = Xmg :: ENDL;
  const Xmg :: t_operation Xmg :: comment    = Xmg :: COMMENT;
  const Xmg :: t_operation Xmg :: clear      = Xmg :: CLEAR;
  const Xmg :: t_operation Xmg :: flush      = Xmg :: FLUSH;
  const Xmg :: t_operation Xmg :: indent     = Xmg :: INDENT;
  const Xmg :: t_operation Xmg :: unindent   = Xmg :: UNINDENT;
  const Xmg :: t_operation Xmg :: addtolast  = Xmg :: ADDTOLAST;
  const Xmg :: t_operation Xmg :: removelast = Xmg :: REMOVELAST;
  const Xmg :: t_operation Xmg :: clearall   = Xmg :: CLEARALL;

  void Xmg :: init (const std::string &_f)
  { 
#ifdef _MPI
    if ( not mpi::main.is_root_node() )
      return;
#endif 
    filename = reformat_home( _f );
    file.open( filename.c_str(), std::ios_base::out|std::ios_base::trunc ); 
    is_empty = true;
    if (file.fail() )
      std::cerr << "Could not open " << filename << std::endl;
    close();
  }
  bool Xmg :: open ()
  {
#ifdef _MPI
    if ( not mpi::main.is_root_node() )
      return false;
#endif
    if ( file.is_open() )
      return true;
    file.open( filename.c_str(), std::ios_base::out|std::ios_base::app ); 
    return file.is_open();
  }    
  void Xmg :: close ()
  {
#ifdef _MPI
    if ( not mpi::main.is_root_node() )
      return;
#endif
    if ( not file.is_open() )
      return;
    flushall();
    file.flush();
    file.close();
  }    
  void Xmg :: flushall ()
  {
#ifdef _MPI
    if ( not mpi::main.is_root_node() )
      return;
#endif 
    if ( line_list.empty() )
      return; 

    open();

    if ( is_empty )  // print revision
      file << "# Subversion Revision " << SVN::Revision << std::endl;

    is_empty = false;

    std::list< std::string > :: const_iterator i_str = line_list.begin();
    std::list< std::string > :: const_iterator i_str_end = line_list.end();
    for(; i_str != i_str_end; ++i_str) 
      file << *i_str << std::endl; 

    line_list.clear();
    file.flush();
  }    

  void Xmg :: add_comment( const std::string &_str )
  {
#ifdef _MPI
    if ( not mpi::main.is_root_node() )
      return;
#endif 
    std::ostringstream sstr;
    sstr << "# ";
    for( types::t_unsigned i = 0; i < indentation; ++i)
      sstr << "  ";
    sstr << _str;
    std::string str = sstr.str();
    line_list.push_back( str );
  }
  void Xmg :: add_to_last( const std::string &_str )
  {
#ifdef _MPI
    if ( not mpi::main.is_root_node() )
      return;
#endif 
    if( line_list.empty() )
      return;
    std::string &str = line_list.back();
    str += _str;
  }
  void Xmg :: add_to_last()
  {
    if( line_list.empty() ) return;

    stream.str(""); stream << line_list.back();
    line_list.pop_back();
  }

  void Xmg :: special_op( t_operation _op )
  {
    switch ( _op )
    {
      case ENDL: line_list.push_back( stream.str() ); stream.str(""); break;
      case CLEAR: stream.str(""); break;
      case COMMENT: stream.str(""); stream << "# "; do_indent(); break;
      case FLUSH: flushall(); break;
      case INDENT: ++indentation; break;
      case UNINDENT: --indentation; break;
      case ADDTOLAST: add_to_last(); break;
      case REMOVELAST: line_list.pop_back(); break;
      case CLEARALL: stream.str(""); clear_all(); break;
    }
  }


  Xmg xmg("convex_hull.agr");
}
