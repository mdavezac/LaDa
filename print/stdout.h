//
//  Version: $Id$
//
#ifndef __PRINT_STDOUT_H_
#define __PRINT_STDOUT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <fstream>
#include <list>

#include "opt/types.h"
#include <revision.h>
#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

#include "operations.h"

namespace Print
{
  //! \brief Simply writes to file. 
  //! \details In mpi codes and production compilation, only the root processor
  //!          can write. In mpi codes with debug enables and the
  //!          _PRINT_ALL_PROCS macro defined, each processor to its own file.
  //!          The first line of the file is always the highest revision number
  //!          found in LaDa. See max_revision.pl in the root directory of
  //!          LaDa.
  //!          StdOut accepts the standard interface of a stream.
  //!  \code
  //     Print::out << "whatever" << some_int << " and " << some_real << Print::endl;
  //!  \endcode
  //!          Note however that the standard formating operations have been
  //!          replaced with Print's very own. See Print::t_Operation for more
  //!          details.
  class StdOut
  {
    protected:
      //! True if nothing has yet been written to file.
      bool is_empty;
      //! True if this processor should write.
      bool do_print;
      //! The name of the file to which to print.
      std::string filename;
      //! The file stream
      std::ofstream file;

    public:
      //! Constructor
      StdOut () : is_empty(true), do_print(false), filename("") {}
      //! Destructor
      ~StdOut() { close(); }
      
      //! Opens the file stream.
      bool open();
      //! Closes the file stream.
      void close();

      //! Returns the filename.
      std::string get_filename() const { return filename; }
      //! Returns true if the filename is not empty.
      bool is_set() const { return not filename.empty(); }
      //! Returns true if the file stream is opened.
      bool is_open() { return do_print and file.is_open(); }
      //! Sets the filename, checks that the file can be opened, and truncates the file.
      void init(const std::string &_f) { if ( filename == _f ) return; init_(_f); }
      //! Handy printing operator in the fashion of standard library streams.
      template<class T_TYPE> inline StdOut& operator<< ( const T_TYPE &_whatever )
        { operator_( _whatever ); return *this; }
      //! A specializable operator<<(). Bit of a hack.
      template<class T_TYPE> inline void operator_ ( const T_TYPE &_whatever );
#ifdef _MPI
#ifndef _PRINT_ALL_PROCS
      //! \ingroup MPI
      //! Syncs filname across all procs.
      void sync_filename() {}
      //! \ingroup MPI
      //! Syncs filname across all procs.
      void sync_filename( std::string & ) {}
#else
      //! \ingroup MPI
      //! Syncs filname across all procs.
      void sync_filename();
      //! \ingroup MPI
      //! Syncs filname across all procs.
      void sync_filename( std::string &_filename );
#endif
#endif

    private:
      //! Initializes new file.
      void init_(const std::string &_f);
      //! \brief Checks that the stream is open.
      //! \details If the file is empty, prints the revision number.
      void do_checks()
      {
        if ( not is_open() ) open();
        if ( not is_empty ) return;
        
        file << "### " << std::endl 
             << "### Subversion Revision Number " << SVN::Revision << std::endl
             << "### " << std::endl;
        is_empty = false;
      }
        
  };

  template<class T_TYPE> inline void StdOut::operator_ ( const T_TYPE &_whatever )
  {
    if ( not do_print ) return;
    do_checks();
    file << _whatever;
  }
  //! Formatting operations with no arguments.
  template<> inline void StdOut::operator_<Print::t_Operation> ( const Print::t_Operation &_op )
  {
    if ( not do_print ) return;
    do_checks();
    Print::apply_ops( file, _op );
  }
  //! Formatting operation with an argument.
  template<> inline void StdOut::operator_<Print::setw> ( const Print::setw &_w )
  {
    if ( not do_print ) return;
    do_checks();
    _w(file);
  }
  //! Formatting operation with an argument.
  template<> inline void StdOut::operator_<Print::setfill> ( const Print::setfill &_w )
  {
    if ( not do_print ) return;
    do_checks();
    _w(file);
  }
  //! Formatting operation with an argument.
  template<> inline void StdOut::operator_<Print::setprecision> ( const Print::setprecision &_w )
  {
    if ( not do_print ) return;
    do_checks();
    _w(file);
  }

  //! Global Print::StdOut object for easy access everywhere.
  extern StdOut out;


}

#endif // _PRINT_XMGRACE_H_
