//
//  Version: $Id$
//


#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python.hpp>
#include <iostream>
#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <stdexcept>

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

using namespace boost::python;

#include "atat.py.hpp"

#define _INMODULE_
BOOST_PYTHON_MODULE(atat)
{
  using namespace atat;
# include "atat.py.hpp"
}

#undef _INMODULE_
