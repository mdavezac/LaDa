#ifdef __DOPYTHON
#ifndef _LADA_PYTHON_CLJ_HPP_
#define _LADA_PYTHON_CLJ_HPP_

#include "LaDaConfig.h"

namespace LaDa
{
  namespace Python
  {
    void expose_clj();
  }
} // namespace LaDa

# ifndef FRIEND_EXPOSE_CLJ
#   define FRIEND_EXPOSE_CLJ friend void LaDa::Python::expose_clj();
# endif
#endif
#else
# ifndef FRIEND_EXPOSE_CLJ
#   define FRIEND_EXPOSE_CLJ 
# endif
#endif
