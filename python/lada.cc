//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/python.hpp>

#include "atat.hpp"
#include "lattice.hpp"
#include "atom.hpp"
#include "structure.hpp"
#include "physics.hpp"
#ifdef __DOCE
#  include "ce.hpp"
#endif

BOOST_PYTHON_MODULE(LaDa)
{
  PythonLaDa::expose_physics();
  PythonLaDa::expose_atat();
  PythonLaDa::expose_lattice();
  PythonLaDa::expose_atom();
  PythonLaDa::expose_structure();
# ifdef __DOCE
  PythonLaDa::expose_ce();
# endif
}
