//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#include <unistd.h>
#include<boost/filesystem/path.hpp>
#include<boost/filesystem/operations.hpp>

#include <iostream>

#include <mpi/mpi_object.h>

#include "fortran_dirio.h"

extern "C" 
{
   //! \brief gets correct path from boost::filesystem::current_path().
   //! \details \a _out is set to the current path and \a _status to 0 if the
   //!          path size is smaller than \a _max_size_out. Otherwise, the
   //!          _max_size_out is set to the path-size and _status set to 1.
   void FC_FUNC_( get_current_directory, GET_CURRENT_DIRECTORY )
                ( char* _out, int *_max_size_out, int *_status )
   {
     const std::string result( boost::filesystem::current_path().string() );
     if( result.size() > *_max_size_out )
     {
       *_max_size_out = result.size();
       *_status = 1u;
       return;
     }
     *_status = 0u;
     std::copy( result.begin(), result.end(), _out );
     std::fill( _out + result.size(), _out + *_max_size_out, ' ' );
   }

   //! \brief Changes current working directory.
   //  \details \a _status is 0 on success, 1 if \a _dir does not exist, 2 if
   //           \a _dir is not a directory, 3 for other errors.
   void FC_FUNC_( change_current_directory, CHANGE_CURRENT_DIRECTORY )
                ( const char* _dir, const int *_size_dir, int *_status )
   {
     if( *_size_dir <= 0 ) { *_status = 3; return; }
     std::string in;
     std::copy( _dir, _dir + *_size_dir, std::back_inserter( in ) );
     const boost::filesystem::path path( in );
     if( not boost::filesystem::exists( path ) )
     {
       *_status = 1u;
       return;
     }
     if( not boost::filesystem::is_directory( path ) )
     {
       *_status = 2u;
       return;
     } 
     chdir( path.string().c_str() );
   }
   
     void FC_FUNC( mpierror, MPIERROR )( int *_ierror )
     {
#      ifdef _MPI
         char error[150];
         int len;
         MPI_Error_string( *_ierror, error, &len );
         std::cerr << error << "\n";
#      endif
     }
}

