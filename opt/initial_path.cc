//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "gsl_lsq.h"

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
  }
} // namespace LaDa
