//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "xml.hpp"

#include <crystal/lattice.h>

namespace PythonLaDa
{
  namespace XML
  {
    template<> std::string nodename<Crystal::Lattice>();
    template<> void do_specialcode< Crystal::Lattice >( Crystal::Lattice &_type );
  }

  void export_lattice();
}
