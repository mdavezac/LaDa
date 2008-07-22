//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/python.hpp>

#include <crystal/atom.h>

namespace PythonLaDa
{
  Crystal::Structure::t_Atom* AtomFromObject( boost::python::list& _ob );
  Crystal::Lattice::t_Site* SiteFromObject( boost::python::list& _ob );

  types::t_real toReal(std::string _str );
  std::string toType( types::t_real _r );
  void export_atom();
}


