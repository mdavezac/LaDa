//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/python.hpp>

#include <revision.h>

#include "atat.hpp"
#include "lattice.hpp"
#include "atom.hpp"
#include "structure.hpp"
#include "physics.hpp"
#ifdef __DOCE
#  include "ce.hpp"
#endif

// namespace PythonLaDa
// {
//   void expose_svn()
//   {
//     using namespace boost::python;
//     def( "Revision", &SVN::Revision );
//     def( "Year",     &SVN::Year );
//     def( "Month",    &SVN::Month );
//     def( "Day",      &SVN::Day );
//     def( "User",     &SVN::User );
//   }
// }

BOOST_PYTHON_MODULE(LaDa)
{
  // PythonLaDa::expose_svn();
  PythonLaDa::expose_physics();
  PythonLaDa::expose_atat();
  PythonLaDa::expose_lattice();
  PythonLaDa::expose_atom();
  PythonLaDa::expose_structure();
# ifdef __DOCE
  PythonLaDa::expose_ce();
# endif
}
