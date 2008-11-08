//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#include <boost/python.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <physics/physics.h>

#include "misc.hpp"
#include "xml.hpp"
#include "physics.hpp"

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
    }

  }
} // namespace LaDa
