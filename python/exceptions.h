#ifndef LADA_PYTHON_EXCEPTIONS_H
#define LADA_PYTHON_EXCEPTIONS_H
#include "LaDaConfig.h"

#include <boost/python/object.hpp>

namespace LaDa
{
  namespace python
  {
    //! \brief Return pointer to root exception class.
    //! \details There is no need to INCREF/DECREF this object. It is
    //!          declared using a static variale and will exist for the
    //!          duration of the program.
    boost::python::object root_exception();
  }
}
# endif 
