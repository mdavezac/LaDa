#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/tuple.hpp>
#include <boost/python/class.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <pyublas/numpy.hpp>

#include <python/std_vector.hpp>

#include "../bandgap.h"
#include "../dipole_elements.h"
#include "dipole_elements.hpp"

namespace LaDa
{
  namespace Python
  {
    namespace bp=boost::python;

    types::t_real oscillator_strength(Pescan::BandGap const &_bg, types::t_real _d, bool _v)
      { return Pescan::oscillator_strength(_bg, _d, _v); }
    types::t_real oscillator_strength2(std::vector<Pescan::Dipole> const &_dipoles)
      { return Pescan::oscillator_strength(_dipoles); }

    pyublas::numpy_vector<types::t_complex> get_p(Pescan::Dipole const &_dipole)
    {
      pyublas::numpy_vector<types::t_complex> result(3);
      result(0) = _dipole.p[0];
      result(1) = _dipole.p[1];
      result(2) = _dipole.p[2];
      return result;
    }
    bp::tuple get_band2band(Pescan::Dipole const &_dipole)
      { return bp::make_tuple(_dipole.band2band.first, _dipole.band2band.second); }

    boost::shared_ptr< std::vector<Pescan::Dipole> >
      get_dipoles(Pescan::BandGap const& _bg, types::t_real _deg) 
      {
        boost::shared_ptr< std::vector<Pescan::Dipole> > result(new std::vector<Pescan::Dipole>);  
        Pescan::dipole_elements(*result, _bg, _deg);
        return result;
      }

    void expose_oscillator_strength()
    {
      namespace bp = boost::python;
      bp::def
      (
        "oscillator_strength",
        &oscillator_strength,
        (
          bp::arg("bandgap"),
          bp::arg("degeneracy") = 0.001,
          bp::arg("verbose") = false
        ),
        "Returns squared norm of the oscillator strength: |<p>|^2.\n\n"
        "  - bandgap is a pescan bandgap functional.\n"
        "  - degeneracy: if |epsilon_i - epsilon_j| < degeneracy, they are"
            "considered a single band. (default: 0.001)\n"
        "  - verbose means verbose output. (default: False)\n" 
      );
      bp::def
      (
        "oscillator_strength",
        &oscillator_strength2,
        bp::arg("dipoles"),
        "Returns squared norm of the oscillator strength: |<p>|^2.\n\n"
        "  - dipoles is a list of dipole elements.\n"
      );
      bp::def 
      (
        "dipole_elements",
        &get_dipoles,
        (
          bp::arg("bangap"),
          bp::arg("degeneracy")
        ),
        "Returns a list of dipole elements from the valence to the conduction band.\n\n"
        "  - bandgap is a pescan bandgap functional.\n"
        "  - degeneracy: if |epsilon_i - epsilon_j| < degeneracy, they are"
            "considered a single band. (default: 0.001)\n"
      );
    }

    void expose_dipoles()
    {
      bp::scope scope = bp::class_<Pescan::Dipole>("Dipole", bp::no_init)
        .def(bp::init<Pescan::Dipole const&>())
        .add_property("p", &get_p)
        .add_property("bands", &get_band2band);

      expose_vector<Pescan::Dipole>("DipoleVector", "List of dipole elements.");
      bp::register_ptr_to_python< boost::shared_ptr< std::vector<Pescan::Dipole> > >();
    }

  } // namespace Python
} // namespace LaDa
