//
//  Version: $Id$
//

#include <boost/python.hpp>
#include <iostream>
#include <algorithm>
#include <boost/lambda/lambda.hpp>

#include <opt/types.h>

#include "atom.h"

#include "atom.py.hpp"

#define _INMODULE_
BOOST_PYTHON_MODULE(atom)
{
   using namespace boost::python;

#  include "atom.py.hpp"
// #  define _TYPE_ types::t_real
// #  define _PYTHONNAME_(object) "r"+ #object
// #  include "atat.py.hpp"
}

#undef _INMODULE_
