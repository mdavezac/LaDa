//
//  Version: $Id$
//
#ifndef __PRINT_XMG_H_
#define __PRINT_XMG_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <fstream>
#include <list>

#include "opt/types.h"
#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

#include "operations.h"

namespace Print
{
  //! Adds comment character to the beginning of each line of a string.
  std::string make_commented_string( const std::string& );



  //! \brief Line based ouputing.
  //! \details This class makes it easy to write one line, go back and modify
  //!          it, delete it... The writing unit is the line.
  //!          In mpi codes and production compilation, only the root processor
  //!          can write. In mpi codes with debug enables and the
  //!          _PRINT_ALL_PROCS macro defined, each processor to its own file.
  //!          The first line of the file is always the highest revision number
  //!          found in LaDa. See max_revision.pl in the root directory of
  //!          LaDa. 
  //!          Xmg accepts the standard interface of a stream.
  //!  \code
  //     Print::xmg << Print::Xmg::comment << "whatever" 
  //                << some_int << " and " << some_real << Print::endl;
  //!  \endcode
  //!          Note however that the standard formating operations have been
  //!          replaced with Print's very own. See Print::t_Operation for more
  //!          details.
  //!          Is also available the following formatting operations:
  //!            - Print::endl goes to the line but does not flush.
  //!            - Print::flush writes to file. It's final.
  //!            - Print::Xmg::comment adds a comment character to the
  //!                                  beginning of the line.
  //!            - Print::Xmg::clearall clears all output which has not yet been
  //!                                   finallized with Print::flush.
  //!            - Print::Xmg::indent adds one more indentation to the current indentation.
  //!            - Print::Xmg::unindent removes one indentation from the current indentation.
  //!            - Print::Xmg::addtolast Goes back to last line (eg undoes all since
  //!                                    last Print::endl call).
  //!            - Print::Xmg::removelast Removes all since penultimuate Print::endl
  //!                                     call  (not including penultimate
  //!                                     Print::endl call).
  //!            - Print::Xmg::clear clears current line.
  //!            .
  class Xmg
  {
    friend std::string make_commented_string( const std::string &_str );

    protected:
      //! Special Print::Xmg operations
      enum t_operation
      { 
        COMMENT,    //!< Adds a comment character
        CLEAR,      //!< Clears current line
        INDENT,     //!< Adds one indentation to next line
        UNINDENT,   //!< Removes one indentation from next line
        ADDTOLAST,  //!< Goes back one line
        REMOVELAST, //!< Removes last line
        CLEARALL    //!< Clears all non-flushed lines.
      };
    public:
      //! Stands for a Print::Xmg standard formatting operation.
      const static t_operation comment;
      //! Stands for a Print::Xmg standard formatting operation.
      const static t_operation clear; 
      //! Stands for a Print::Xmg standard formatting operation.
      const static t_operation indent;
      //! Stands for a Print::Xmg standard formatting operation.
      const static t_operation unindent; 
      //! Stands for a Print::Xmg standard formatting operation.
      const static t_operation addtolast;
      //! Stands for a Print::Xmg standard formatting operation.
      const static t_operation removelast;
      //! Stands for a Print::Xmg standard formatting operation.
      const static t_operation clearall;
    protected:
      //! The comment charracter
      const static std::string comment_string;

    protected:
      //! The number of current indentation
      types::t_unsigned indentation;
      //! The name of the file to which to write
      std::string filename;
      //! A stream for the file to which to write
      std::ofstream file;
      //! \brief The lines which have been written (Print::endl) but no flushed
      //!        to file.
      std::list< std::string > line_list;
      //! A stream for writting the current line
      std::ostringstream stream;
      //! True if nothing has yet been written to file.
      bool is_empty;
      //! Wether this processor should write to file.
      bool do_print;

    public:
      //! Constructor
      Xmg() : indentation(0), filename(""), is_empty(true), do_print(false) {};
      //! Destructor
      ~Xmg() { close(); }
      
      //! Opens the file and stream.
      bool open();
      //! Closes the file and stream.
      void close();

      //! Clears all non-flushed lines.
      void clear_all()
      {
        if ( not do_print )  return;
        line_list.clear(); 
      }
      //! Opens the file stream.
      bool is_open() 
      {
        return do_print and file.is_open(); 
      }
      //! Adds a line to the list of lines. Does not flush to file.
      void add_line( const std::string &_str )
      {
        if ( not do_print )  return;
        line_list.push_back(_str); 
      }
      //! Adds a commented line to the list of lines. Does not flush to file.
      void add_comment( const std::string &_str);
      //! Appends \a _str to the last line.
      void add_to_last( const std::string &_str);
      //! Removes last line.
      void remove_last()
      {
        if ( not do_print )  return;
        line_list.pop_back(); 
      }
      //! Initializes the file.
      void init(const std::string &_f);
      //! Returns the filename
      std::string get_filename() const { return filename; }
      //! Easy printing in standard fashion.
      template<class T_TYPE> Xmg& operator<< (const T_TYPE &_whatever);

    protected:
      //! Flushes all lines to file and flushed file.
      void flushall();
      //! Writes to current line.
      template<class T_TYPE> void to_current_stream ( const T_TYPE& _type )
        { if( not do_print ) return; stream << _type; }
      //! Adds indentation to current line.
      void do_indent()
        { for( types::t_unsigned i = 0; i < indentation; ++i) stream << "  "; }
      //! Removes current line and goes back to last.
      void add_to_last();
      //! Specializable printing operation.
      void special_op( t_operation _op);


#ifdef _MPI
    public:
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
  };

  //! \brief Print::Xmg object for easy access everywhere.
  extern Xmg xmg;

  //! Specialized Print::Xmg formatting operations for current line.
  template<> inline void Xmg :: to_current_stream<Xmg::t_operation>(const Xmg::t_operation &_op)
  {
    special_op( _op );
  }
  //!  Specialized standard formatting operations for current line.
  template<> inline void Xmg :: to_current_stream<Print::t_Operation>(const Print::t_Operation &_op)
  {
    if ( _op == Print::ENDL  and stream.str().empty() ) return;
    switch ( _op )
    {
      case ENDL: line_list.push_back( stream.str() ); stream.str(""); break;
      case FLUSH: flushall(); break;
      default: Print::apply_ops( stream, _op ); break;
    }
  }
  //!  Specialized standard formatting operations for current line.
  template<> inline void Xmg::to_current_stream<Print::setw>(const Print::setw &_op) { _op(stream); }
  //!  Specialized standard formatting operations for current line.
  template<> inline void Xmg::to_current_stream<Print::setfill>(const Print::setfill &_op) { _op(stream); }
  //!  Specialized standard formatting operations for current line.
  template<> inline void Xmg::to_current_stream<Print::setprecision>(const Print::setprecision &_op)
    { _op(stream); }
   

  template<class T_TYPE> inline Xmg& Xmg::operator<< (const T_TYPE &_whatever)
  {
    if ( do_print ) to_current_stream( _whatever );
    return *this;
  }

}

#endif // _PRINT_XMGRACE_H_
