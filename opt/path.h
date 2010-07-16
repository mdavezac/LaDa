#ifndef LADA_PATH_H
#define LADA_PATH_H

#include "LaDaConfig.h"

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
