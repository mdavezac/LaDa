#ifndef LADA_STDOUT_H
#define LADA_STDOUT_H
#include "LaDaConfig.h"
#include "FCMangle.h"

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <fstream>
#include <boost/utility.hpp>
#include <boost/filesystem/path.hpp>

// declares fortran interface
//! \cond
extern "C"
{
  void FC_GLOBAL_(lada_redirect_open,  LADA_REDIRECT_OPEN )( const int *, const int *,
                                                             const char *, int *, const int *);
# ifdef FUCKING_CRAY
    void FC_GLOBAL_(fucking_cray, FUCKING_CRAY)(int *);
# else
    void FC_GLOBAL_(lada_redirect_close, LADA_REDIRECT_CLOSE)(int *);
# endif
}
//! \endcond

namespace LaDa
{
  namespace opt
  {
    //! \brief Redirects file descriptor to other file. 
    //! \details In practice this class can be used to redirect the standard output
    //! and standard error to actual files. The redirection starts at
    //! construction time and stops after destruction, eg from declaration to
    //! out-of-scope.
    class Redirect : boost::noncopyable
    {
      public:
        //! Construction
        Redirect(int const _which, boost::filesystem::path const &_out, bool _append = false)
        {
          which_ = _which;
          if(not _append)
            filedesc_ = ::open(_out.string().c_str(), O_CREAT|O_TRUNC |O_WRONLY, S_IWUSR|S_IRUSR);
          else           
            filedesc_ = ::open(_out.string().c_str(), O_CREAT|O_APPEND|O_WRONLY, S_IWUSR|S_IRUSR);
          if(filedesc_  < 0) return; 
          save_  = dup(_which);
  	  dup2(filedesc_, which_);
        }
        //! Destruction;
        virtual ~Redirect() { close(); }

        //! Closes the pipe.
        void close() 
        {
          if( filedesc_ < 0 ) return;
          ::fsync(filedesc_);
          ::close(filedesc_);
          filedesc_ = -1;
  	  dup2(save_, which_);
        }
        //! True if a pipe is open.
        bool is_open() const { return filedesc_ >= 0; }

      private:
        //! Which file descriptor to pipe from.
        int which_;
        //! Which file descriptor to pipe to.
        int save_;
        //! The file descriptor where to save.
        int filedesc_;
    };

    //! \brief Redirects fortran standard output to file.
    //! \details  The redirection starts at construction time and stops after
    //! destruction, eg from declaration to out-of-scope. If the path is empty,
    //! then redirects to /dev/null.
    class FortranRedirect : boost::noncopyable
    {
      public:
        //! Which unit to redirect.
        enum Unit 
        {
          input  = 5, //! redirects input unit.
          output = 6, //! redirects input unit.
          error  = 0  //! redirects input unit.
        };
        //! Construction
        FortranRedirect( Unit const _which, boost::filesystem::path const &_out,
                         bool _append = false)
        {
          std::string const path = _out.string().size() ? _out.string(): "/dev/null";
          int const N = path.size();
          which_ = _which;
          int isopen = 1;
          int const append = _append ? 1: 0;
          FC_GLOBAL_(lada_redirect_open, LADA_REDIRECT_OPEN)
            (&which_, &N, path.c_str(), &isopen, &append);
          isopen_ = isopen == 1;
        }
        //! Destruction;
        virtual ~FortranRedirect() { close(); }

        //! Closes the pipe.
        void close() 
        {
          if( not isopen_ ) return;
          isopen_ = false;
#         ifdef FUCKING_CRAY
            FC_GLOBAL_(fucking_cray, FUCKING_CRAY)(&which_);
#         else
            FC_GLOBAL_(lada_redirect_close, LADA_REDIRECT_CLISE)(&which_);
#         endif
        }
        //! True if a pipe is open.
        bool is_open() const { return isopen_; }

      private:
        //! Which file descriptor to pipe from.
        int which_;
        //! The file descriptor where to save.
        bool isopen_;
    };
  }
}

#endif 
