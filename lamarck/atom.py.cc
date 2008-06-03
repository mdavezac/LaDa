//
//  Version: $Id$
//

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <iostream>
#include <algorithm>
#include <boost/lambda/lambda.hpp>

#include <opt/types.h>
#include <opt/convex_hull.h>
#include <opt/debug.h>

#include "functional_builder.h"
#include "constituent_strain.h"
#include "harmonic.h"
#include "lattice.h"
#include "structure.h"
#include "atom.h"

#include<tinyxml/tinyxml.h>

using namespace boost::python;
#include "atom.py.hpp"

#define _INMODULE_
BOOST_PYTHON_MODULE(atom)
{

#  include "atom.py.hpp"
// #  define _TYPE_ types::t_real
// #  define _PYTHONNAME_(object) "r"+ #object
// #  include "atat.py.hpp"
}

#undef _INMODULE_
