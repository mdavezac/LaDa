//
//  Version: $Id$
//
#ifndef _OPT_INITIAL_PATH_H_
#define _OPT_INITIAL_PATH_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/filesystem/path.hpp>

#include "debug.h"

namespace LaDa
{
  namespace opt
  {
    //! Keeps path when entering main().
    class InitialPath
    {
      public:
        //! Initializes the path.
        static void init();
        //! Return a reference to the path.
        static boost::filesystem::path& path()
        { 
          __ASSERT( is_initialized_ == false, "Path has not been initialized.\n" )
          return path_; 
        }

      protected:
        //! The path on call to InitialPath.
        static boost::filesystem::path path_;
        //! Whether InitialPath::init() has been called.
        static bool is_initialized_;
    };

  }
} // namespace LaDa

#endif
