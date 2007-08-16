#include<iostream>
#include <sstream>

#include "print_xmgrace.h"
#include "pescan_interface/interface.h"
namespace darwin
{

  void PrintXmg :: init (const std::string &_f)
  { 
#ifdef _MPI
    if ( not mpi::main.is_root_node() )
      return;
#endif 
    filename = reformat_home( _f );
    file.open( filename.c_str(), std::ios_base::out|std::ios_base::trunc ); 
    if (file.fail() )
      std::cerr << "Could not open " << filename << std::endl;
    close();
  }
  bool PrintXmg :: open ()
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
  void PrintXmg :: close ()
  {
#ifdef _MPI
    if ( not mpi::main.is_root_node() )
      return;
#endif
    if ( not file.is_open() )
      return;
    flush();
    file.flush();
    file.close();
  }    
  void PrintXmg :: flush ()
  {
#ifdef _MPI
    if ( not mpi::main.is_root_node() )
      return;
#endif 
    if ( line_list.empty() )
      return; 

    open();

    std::list< std::string > :: const_iterator i_str = line_list.begin();
    std::list< std::string > :: const_iterator i_str_end = line_list.end();
    for(; i_str != i_str_end; ++i_str) 
      file << *i_str << std::endl; 

    line_list.clear();
    file.flush();
  }    

  void PrintXmg :: add_comment( const std::string &_str )
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
  void PrintXmg :: add_to_last( const std::string &_str )
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


  PrintXmg printxmg("convex_hull.agr");
}
