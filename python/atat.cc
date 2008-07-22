//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python.hpp>
#include <iostream>
#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <stdexcept>

#include <opt/types.h>

#define STREAM_VECTOR
#include <atat/fxvector.h>
#include <atat/vectmac.h>

#ifdef _TYPE_
#error "Please change _TYPE_ to something else, as it is already being used."
#endif
#ifdef _DIM_
#error "Please change _TYPE_ to something else, as it is already being used."
#endif

using namespace boost::python;

namespace PythonLaDa
{
  using namespace atat;
# include "atat.hpp"

# define _INMODULE_
  void export_atat()
  {
    using namespace atat;
#   include "atat.hpp"
  }

#undef _INMODULE_
}
