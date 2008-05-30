//
//  Version: $Id$
//

// #ifdef HAVE_CONFIG_H
// #include <config.h>
// #endif

#include <boost/python.hpp>
#include <iostream>
#include <algorithm>
#include <boost/lambda/lambda.hpp>

// #include <boost/python/module.hpp>
// #include <boost/python/def.hpp>

// 
//  char const* greet()
//  {
//     return "hello, world";
//  }
// 
//  BOOST_PYTHON_MODULE(atat)
//  {
//      using namespace boost::python;
//      def("greet", greet);
//  }
// 
#include <opt/types.h>

#define STREAM_VECTOR
#include "fxvector.h"
#include "vectmac.h"

#ifdef _TYPE_
#error "Please change _TYPE_ to something else, as it is already being used."
#endif
#ifdef _DIM_
#error "Please change _TYPE_ to something else, as it is already being used."
#endif

#  include "atat.py.hpp"

#define _INMODULE_
BOOST_PYTHON_MODULE(atat)
{
  using namespace boost::python;
  using namespace atat;
#  include "atat.py.hpp"
// #  define _TYPE_ types::t_real
// #  define _PYTHONNAME_(object) "r"+ #object
// #  include "atat.py.hpp"
}

#undef _INMODULE_
