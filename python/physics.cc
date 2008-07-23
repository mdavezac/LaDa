//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <boost/python.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <opt/physics.h>

#include "misc.hpp"
#include "xml.hpp"
#include "physics.hpp"

namespace PythonLaDa
{
  void expose_physics()
  {
    using namespace boost::python;
    BOOST_PYTHON_MODULE( Physics )
    {
      def("Z", Physics::Z);
      def("Symbol", Physics::Symbol);
      types::t_unsigned (*ptr_charge) (const std::string &) = &Physics::Charge;
      def("Charge", ptr_charge);
      def("Mass", Physics::Mass);
    }
  }

}
