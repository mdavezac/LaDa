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
      def("Z", Physics::Atomic::Z,
          "Given an atomic symbol, returns the atomic number.");
      def("Symbol", Physics::Atomic::Symbol,
          "Given an atomic number, returns the atomic symbol.");
      types::t_unsigned (*ptr_charge) (const std::string &) = &Physics::Atomic::Charge;
      def("Charge", ptr_charge,
          "Given an atomic symbol, returns the number of valence electrons." );
      def("Mass", Physics::Atomic::Mass,
          "Given an atomic symbol, returns the atomic mass.");
      def("a0", Physics::a0,
          "Returns the Bhor radius in A, nm, m, or cm" );
      def("Hartree", Physics::Hartree,
          "Returns the Hartree energy in eV, Rydberg, or Hartree" );
      def("Rydberg", Physics::Rydberg,
          "Returns the Rydberg energy in eV, Rydberg, or Hartree" );
      def("emass", Physics::emass, "Returns the mass of the electron in eV, amu, kg, g, MeV." );
      def("hbar", Physics::hbar, "Returns hbar in  eV*s, erg*s, J*s." );
      def("planck", Physics::Planck, "Returns the Planck constant in  eV*s, erg*s, J*s, Ry, H." );
    }

  }
} // namespace LaDa

BOOST_PYTHON_MODULE(physics)
{
  LaDa::Python::expose_physics();
}
