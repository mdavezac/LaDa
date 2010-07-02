#ifndef LADA_PATH_H
#define LADA_PATH_H

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/filesystem/path.hpp>

namespace LaDa
{
  namespace opt
  {
    //! Expand all user home paths and environment variables.
    boost::filesystem::path const expand_path(boost::filesystem::path const &_path);
  }
}
#endif
