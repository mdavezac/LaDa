//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/filesystem/operations.hpp>

#include "initial_path.h"


namespace LaDa
{
  namespace opt
  {
    void InitialPath :: init()
    {
      __ASSERT( is_initialized_ == true,
                "opt::InitialPath::init has already been called.\n" )
      path_ =  boost::filesystem::initial_path();
      is_initialized_ = true;
    }

    boost::filesystem::path InitialPath::path_;
    bool InitialPath::is_initialized_=false;

    bool create_directory(boost::filesystem::path const &_path)
    {
      namespace bfs = boost::filesystem;
      if( bfs::exists(_path) ) return bfs::is_directory(_path);
      if( _path.has_parent_path() )
        if(not opt::create_directory(_path.parent_path())) return false;
      // creates directory.
      bfs::create_directory(_path);
      return bfs::exists(_path) and bfs::is_directory(_path);
    }
  }
} // namespace LaDa
