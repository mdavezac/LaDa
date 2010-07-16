#ifdef LADA_DO_PYTHON
#ifndef _LADA_PYTHON_CLJ_HPP_
#define _LADA_PYTHON_CLJ_HPP_

#include "LaDaConfig.h"

namespace LaDa
{
  namespace python
  {
    void expose_clj();
  }
} // namespace LaDa

# ifndef FRIEND_EXPOSE_CLJ
#   define FRIEND_EXPOSE_CLJ friend void LaDa::python::expose_clj();
# endif
#endif
#else
# ifndef FRIEND_EXPOSE_CLJ
#   define FRIEND_EXPOSE_CLJ 
# endif
#endif
