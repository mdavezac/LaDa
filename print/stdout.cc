//
//  Version: $Id$
//
#include<iostream>
#include <sstream>

#include <revision.h>

#include "stdout.h"
#include "manip.h"

namespace Print
{

  void StdOut :: init_ (const std::string &_f)
  { 
    if ( is_open() ) close();
    filename = reformat_home( _f );
    do_print = true;
#ifdef _MPI
    if ( not mpi::main.is_root_node() )
    {
      std::ostringstream sstr; 
      sstr << filename << ".mpi:" <<  mpi::main.rank();
      filename = sstr.str();
    }
#ifndef _PRINT_ALL_PROCS
    do_print = mpi::main.is_root_node();
#endif
#endif 
    is_empty = true;
    if ( not do_print ) return;
    file.open( filename.c_str(), std::ios_base::out|std::ios_base::trunc ); 
    if (file.fail() )
    {
      std::cerr << "Could not open " << filename << std::endl;
      do_print = false;
      return;
    }
    close();
  }
  bool StdOut :: open ()
  {
    if ( not do_print ) return true;
    if ( file.is_open() ) return true;
    file.open( filename.c_str(), std::ios_base::out|std::ios_base::app ); 
    return file.is_open();
  }    
  void StdOut :: close ()
  {
    if ( not do_print ) return;
    if ( not file.is_open() ) return;
    file.flush();
    file.close();
  }    

#ifdef _MPI
#ifdef _PRINT_ALL_PROCS
  void StdOut::sync_filename( std::string &_filename )
  {
    filename = _filename;
    sync_filename();
  }
  void StdOut::sync_filename()
  {
    mpi::BroadCast bc( mpi::main );
    bc << filename
       << mpi::BroadCast::allocate
       << filename
       << mpi::BroadCast::broadcast
       << filename;

    init_(filename);
  }
#endif
#endif

  StdOut out;
}
