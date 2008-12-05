//
//  Version: $Id$
//
#ifndef _LADA_FORTRAN_DIRIO_H_
#define _LADA_FORTRAN_DIRIO_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

extern "C" 
{
   //! \brief gets correct path from boost::filesystem::current_path().
   //! \details \a _out is set to the current path and \a _status to 1 if the
   //!          path size is smaller than \a _max_size_out. Otherwise, the
   //!          _max_size_out is set to the path-size and _status set to zero.
   void FC_FUNC_( get_current_directory, GET_CURRENT_DIRECTORY )
                ( char* _out, int *_max_size_out, int *_status );
   //! \brief Changes current working directory.
   //  \details \a _status is 0 on success, 1 if \a _dir does not exist, 2 if
   //           \a _dir is not a directory.
   void FC_FUNC_( change_current_directory, CHANGE_CURRENT_DIRECTORY )
                ( const char* _dir, const int *_max_size_out, int *_status );
}
#endif
