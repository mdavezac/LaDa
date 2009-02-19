//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#include <boost/python.hpp>

#include <opt/types.h>
#include <opt/debug.h>

#include <python/misc.hpp>
#include <python/xml.hpp>

#include "../physics.h"

namespace LaDa
{
  namespace Python
  {
    void expose_physics()
    {
      using namespace boost::python;
      def("Z", Physics::Atomic::Z);
      def("Symbol", Physics::Atomic::Symbol);
      types::t_unsigned (*ptr_charge) (const std::string &) = &Physics::Atomic::Charge;
      def("Charge", ptr_charge);
      def("Mass", Physics::Atomic::Mass);
      def("a0", Physics::a0);
      def("Hartree", Physics::Hartree);
      def("Rydberg", Physics::Rydberg);
    }

  }
} // namespace LaDa

BOOST_PYTHON_MODULE(Physics)
{
  LaDa::Python::expose_physics();
}
