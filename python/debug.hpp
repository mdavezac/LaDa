#ifndef _LADA_PYTHON_DEBUG_HPP_
#define _LADA_PYTHON_DEBUG_HPP_

#include "LaDaConfig.h"

#include <sstream>
#include <string>
#include <boost/python/errors.hpp>

#define LADA_PYTHON_ERROR( error_code, error_string ) \
  {\
    std::ostringstream sstr; \
    sstr << (std::string(__FILE__) + ", line: ") <<  __LINE__ << ":\n" \
         << error_string;\
    std::cerr << sstr.str() << "\n"; \
    PyErr_SetString( error_code, sstr.str().c_str() );\
  }

#endif // _LADA_PYTHON_DEBUG_HPP_
