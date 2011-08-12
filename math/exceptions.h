#ifndef LADA_LNS_EXCEPTIONS_H_
#define LADA_LNS_EXCEPTIONS_H_

#include "LaDaConfig.h"

#include <root_exceptions.h>

namespace LaDa
{
  namespace error
  {
    //! Root of input errors.
    struct array_of_different_sizes: virtual internal {};
  }
}

#endif
