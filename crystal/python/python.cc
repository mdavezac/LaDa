//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>

#include "structure.hpp"
#include "atom.hpp"

BOOST_PYTHON_MODULE(Crystal)
{
  LaDa::Python::expose_atom();
  LaDa::Python::expose_structure();
}
