#include "print_xmgrace.h"

#include<iostream>
#include <sstream>

#include<mpi/mpi_object.h>

namespace darwin
{

  void PrintXmg :: init (const std::string &_f)
  { 
    filename = _f;
    if ( mpi::main.size() > 1 )
    {
      std::ostringstream sstr;
      sstr << filename << "." << mpi::main.rank();
      filename = sstr.str();
    }
    file.open( filename.c_str(), std::ios_base::out|std::ios_base::trunc ); 
    close();
  }
  bool PrintXmg :: open ()
  {
    if ( file.is_open() )
      return true;;
    file.open( filename.c_str(), std::ios_base::out|std::ios_base::app ); 
    return file.is_open();
  }    
  void PrintXmg :: close ()
  {
    if ( not file.is_open() )
      return;
    flush();
    file.flush();
    file.close();
  }    
  void PrintXmg :: flush ()
  {
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
    if( line_list.empty() )
      return;
    std::string &str = line_list.back();
    str += _str;
  }


  PrintXmg printxmg("convex_hull.agr");
}
