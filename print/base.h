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
#include <boost/filesystem/path.hpp>

#include <opt/types.h>

#include <mpi/mpi_object.h>

namespace boost {
  namespace serialization {

    //! Serializes a path.
    template<class Archive>
    void serialize(Archive & _ar, boost::filesystem::path & _p, const unsigned int version)
    {
      std::string path( _p.string() );
      _ar & path; 
      _p = path;
    }
  }
}
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
    public:
      //! Type of the path.
      typedef boost::filesystem::path t_Path;
    protected:
      //! True if nothing has yet been written to file.
      bool is_empty;
      //! True if this processor should write.
      bool do_print;
      //! Whether to truncate or append;
      bool truncate;
      //! The name of the file to which to print.
      t_Path filename;
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
      const t_Path& get_filename() const { return filename; }
      //! Returns true if the filename is not empty.
      bool is_set() const { return not filename.empty(); }
      //! Returns true if the file stream is opened.
      bool is_open() { return do_print and file.is_open(); }
      //! Sets the filename, checks that the file can be opened, and truncates the file.
      void init(const t_Path &_f);
#ifdef _MPI
      //! \ingroup MPI
      //! Syncs filname across all procs.
      void sync_filename();
      //! \ingroup MPI
      //! Syncs filname across all procs.
      void sync_filename( t_Path &_filename );
#endif
      //! Sets to appending file when running init.
      void dont_truncate() { truncate = false; }
      //! Sets to truncating file when running init.
      void do_truncate() { truncate = true; }
      //! Sets to printing.
      void doprint( bool _d  ) { do_print = _d; init_( filename ); }

    private:
      //! Initializes new file.
      void init_(const t_Path &_f);
      //! \brief Checks that the stream is open.
      //! \details If the file is empty, prints the revision number.
      void do_checks();
  };

}

#endif // _PRINT_XMGRACE_H_
