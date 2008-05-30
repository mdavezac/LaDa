//
//  Version: $Id$
//

#include <boost/python.hpp>
#include <iostream>
#include <algorithm>
#include <boost/lambda/lambda.hpp>

#include <opt/types.h>

#include "structure.h"
#include "atom.h"

typedef Ising_CE::Structure::t_Atom t_Atom;

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
