//
//  Version: $Id$
//
#include<iostream>
#include <sstream>

#include <revision.h>

#include "base.h"
#include "manip.h"

namespace Print
{
  void Base :: do_checks()
  {
    if ( not is_open() ) open();
    if ( not is_empty ) return;
    
    file << "### " << std::endl 
         << "### Subversion Revision Number " << SVN::Revision << std::endl
         << "### " << std::endl;
    is_empty = false;
  }

  void Base :: init_ (const std::string &_f)
  { 
    if ( is_open() ) close();
    filename = reformat_home( _f );
    do_print = true;
    __NOTMPIROOT( 
      std::ostringstream sstr; 
      sstr << filename << ".mpi:" <<  mpi::main.rank();
      filename = sstr.str();
    )
    do_print = false;
    __ROOTCODE( do_print = true; )
    is_empty = true;
    if ( not do_print ) return;
    
    file.open( filename.c_str(),
               truncate ? std::ios_base::trunc: std::ios_base::app | std::ios_base::out); 

    if (file.fail() )
    {
      std::cerr << "Could not open " << filename << std::endl;
      do_print = false;
      return;
    }
    if( not truncate ) file << "\n\n"; 
    close();
  }
  bool Base :: open ()
  {
    if ( not do_print ) return true;
    if ( file.is_open() ) return true;
    file.open( filename.c_str(), std::ios_base::out|std::ios_base::app ); 
    return file.is_open();
  }    
  void Base :: close ()
  {
    if ( not do_print ) return;
    if ( not file.is_open() ) return;
    file.flush();
    file.close();
  }    

#ifdef _MPI
  void Base::sync_filename( std::string &_filename )
  {
    filename = _filename;
    sync_filename();
    _filename = filename;
  }
  void Base::sync_filename()
  {
    boost::mpi::broadcast( mpi::main, filename, 0 );
    init_(filename);
  }
#endif

}
