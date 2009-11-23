//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include<iostream>
#include <sstream>
#include<boost/lexical_cast.hpp>
#include<boost/filesystem/operations.hpp>


#include "base.h"
#include "manip.h"

#ifdef _MPI
# include <boost/mpi/collectives.hpp>
#endif

namespace LaDa
{
  namespace Print
  {
    void Base :: do_checks()
    {
      if ( not is_open() ) open();
      if ( not is_empty ) return;
      
      is_empty = false;
    }

    void Base :: init(const t_Path &_f)
    { 
      do_print = false;
      __MPICODE( boost::mpi::communicator world; )
      __ROOTCODE( world, do_print = true; )
      if ( filename == _f ) return;
      init_(_f); 
    }

    void Base :: init_ (const t_Path &_f)
    { 
      if ( is_open() ) close();
      filename = reformat_home( _f.string() );
      do_print = true;
      __MPICODE( boost::mpi::communicator world; )
      __NOTMPIROOT( world, 
        filename =    filename.string()
                    + ".mpi:" 
                    + boost::lexical_cast<std::string>( world.rank() );
      )
      is_empty = true;
      if ( (not do_print) or (not boost::filesystem::exists( filename )) ) return;

      if( truncate ) { boost::filesystem::remove( filename ); return; }

      file.open( filename.string().c_str(),
                 truncate ? std::ios_base::trunc: std::ios_base::app | std::ios_base::out); 
      if (file.fail() )
      {
        std::cerr << "Could not open " << filename << std::endl;
        do_print = false;
        return;
      }
      file << "\n\n"; 
      close();
    }
    bool Base :: open ()
    {
      if ( not do_print ) return true;
      if ( file.is_open() ) return true;
      file.open( filename.string().c_str(), std::ios_base::out|std::ios_base::app ); 
      return file.is_open();
    }    
    void Base :: close ()
    {
      if ( not do_print ) return;
      if ( not file.is_open() ) return;
      file.flush();
      file.close();
    }    

#   ifdef _MPI
      void Base::sync_filename( t_Path &_filename )
      {
        filename = _filename;
        sync_filename();
        _filename = filename;
      }
      void Base::sync_filename()
      {
        boost::mpi::communicator world;
        boost::mpi::broadcast( world, filename, 0 );
        init_(filename);
      }
#   endif

  }
} // namespace LaDa
