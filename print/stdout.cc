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
#ifdef _MPI
    if ( not mpi::main.is_root_node() )
    {
      std::ostringstream sstr; 
      sstr << filename << ".mpi:" <<  mpi::main.rank();
      filename = sstr.str();
    }
#endif 
    is_empty = true;
    if ( not do_print ) return;
    file.open( filename.c_str(), std::ios_base::out|std::ios_base::trunc ); 
    if (file.fail() ) std::cerr << "Could not open " << filename << std::endl;
    close();
  }
  bool StdOut :: open ()
  {
    if ( not do_print ) return false;
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

#if defined(_MPI) and defined(_PRINT_ALL_PROCS)
  void StdOut::sync_filename()
  {
    mpi::BroadCast bc( mpi::main );
    bc << filename
       << mpi::BroadCast::allocate
       << filename
       << mpi::BroadCast::broadcast
       << filename;

    if ( mpi::main.is_root_node() ) return;

    init_(filename);
  }
#endif

  StdOut out("out");
}
