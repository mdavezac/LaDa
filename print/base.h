//
//  Version: $Id$
//
#ifndef __PRINT_BASE_H_
#define __PRINT_BASE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <fstream>

#include <opt/types.h>

#include <mpi/mpi_object.h>

namespace Print
{
  //! \brief Base class for printing to file
  //! \details This is a base class for printing to file. Mainly, it prints the
  //!          revision number at the beginning of the file. It also handles
  //!          mpi. By default, only the root node can print. When compiled
  //!          with the \a _PRINT_ALL_PROCS defined, each node print to its own
  //!          file.
  class Base
  {
    protected:
      //! True if nothing has yet been written to file.
      bool is_empty;
      //! True if this processor should write.
      bool do_print;
      //! Whether to truncate or append;
      bool truncate;
      //! The name of the file to which to print.
      std::string filename;
      //! The file stream
      std::ofstream file;

    public:
      //! Constructor
      Base () : is_empty(true), do_print(false), truncate(true), filename("") {}
      //! Destructor
      ~Base() { close(); }
      
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
#ifdef _MPI
      //! \ingroup MPI
      //! Syncs filname across all procs.
      void sync_filename();
      //! \ingroup MPI
      //! Syncs filname across all procs.
      void sync_filename( std::string &_filename );
#endif
      //! Sets to appending file when running init.
      void dont_truncate() { truncate = false; }
      //! Sets to truncating file when running init.
      void do_truncate() { truncate = true; }

    private:
      //! Initializes new file.
      void init_(const std::string &_f);
      //! \brief Checks that the stream is open.
      //! \details If the file is empty, prints the revision number.
      void do_checks();
      void doprint( bool _d  ) { do_print = _d; init( filename ); }
  };

}

#endif // _PRINT_XMGRACE_H_
