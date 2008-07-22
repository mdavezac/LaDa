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
#include "ce.hpp"

BOOST_PYTHON_MODULE(LaDa)
{
  export_atat();
  export_lattice();
  export_atom();
  export_structure();
# ifdef __DOCE
  export_ce();
# endif
}
