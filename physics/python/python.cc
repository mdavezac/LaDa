#include "LaDaConfig.h"


#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>

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
      namespace bp = boost::python;
      bp::scope scope;
      scope.attr("__doc__") = "Package with text-book physical quantities.";
      bp::docstring_options doc_options(true, false);

      bp::def("Z", Physics::Atomic::Z,
              "Given an atomic symbol, returns the atomic number.");
      bp::def("Symbol", Physics::Atomic::Symbol,
              "Given an atomic number, returns the atomic symbol.");
      types::t_unsigned (*ptr_charge) (const std::string &) = &Physics::Atomic::Charge;
      bp::def("Charge", ptr_charge,
              "Given an atomic symbol, returns the number of valence electrons." );
      bp::def("Mass", Physics::Atomic::Mass,
              "Given an atomic symbol, returns the atomic mass.");
      bp::def("a0", Physics::a0,
              "Returns the Bhor radius in A, nm, m, or cm" );
      bp::def("Hartree", Physics::Hartree,
              "Returns the Hartree energy in eV, Rydberg, or Hartree" );
      bp::def("Rydberg", Physics::Rydberg,
              "Returns the Rydberg energy in eV, Rydberg, or Hartree" );
      bp::def("emass", Physics::emass, "Returns the mass of the electron in eV, amu, kg, g, MeV." );
      bp::def("hbar", Physics::hbar, "Returns hbar in  eV*s, erg*s, J*s." );
      bp::def("planck", Physics::Planck, "Returns the Planck constant in  eV*s, erg*s, J*s, Ry, H." );
    }

  }
} // namespace LaDa

BOOST_PYTHON_MODULE(physics)
{
  LaDa::Python::expose_physics();
}
