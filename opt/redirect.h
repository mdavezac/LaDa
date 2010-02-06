#ifndef LADA_STDOUT_H
#define LADA_STDOUT_H

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <fstream>
#include <boost/utility.hpp>
#include <boost/filesystem/path.hpp>


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
  }
}

#endif 
