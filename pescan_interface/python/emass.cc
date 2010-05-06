//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/tuple.hpp>
#include <boost/python/list.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

#include "../emass.h"
#include "emass.hpp"

namespace LaDa
{
  namespace python 
  {
    struct eMass : public LaDa::Pescan::eMass
    {
      types::t_real tolerance;
      types::t_unsigned itermax;
      bool verbose;

      eMass() : tolerance(1e-4), itermax(10), verbose(false) {}
      eMass   ( const eMass &_c )
            : Pescan::eMass( _c ),
              tolerance( _c.tolerance ),
              itermax( _c.itermax ),
              verbose( _c.verbose ) {}
      math::rVector3d get_direction() const { return direction; }
      void set_direction(math::rVector3d const &_direction) { direction = _direction; }
      math::rVector3d get_kpoint() const { return kpoint; }
      void set_kpoint(math::rVector3d const &_kpoint) { kpoint = _kpoint; }
    };

    
    boost::python::list get_masses
    (
      eMass& _emass,
      const Pescan::Interface& _interface,
      const math::rMatrix3d &_ocell, 
      const Crystal::Structure &_structure,
      const types::t_real &_eref
    )
    {
      typedef std::vector< std::pair< types::t_real, types::t_real > > t_Out;
      t_Out out;
      _emass.cgs.tolerance = _emass.tolerance;
      _emass.cgs.verbose = _emass.verbose;
      _emass.cgs.itermax = _emass.itermax;
      _emass
      (
       _interface,
       _ocell, 
       _structure,
       _eref,
       out 
      );
      boost::python::list result;
      foreach( const t_Out :: value_type &r, out )
        result.append( boost::python::make_tuple( r.first, r.second ) );
      return result;
    }
    boost::python::list get_masses2
    (
      eMass& _emass,
      const Pescan::Interface& _interface,
      const math::rMatrix3d &_ocell, 
      const Crystal::TStructure<std::string> &_structure,
      const types::t_real &_eref
    )
    {
      Crystal::Structure structure;
      Crystal::convert_string_to_real_structure( _structure, structure );     
      return get_masses(_emass, _interface, _ocell, structure, _eref);
    }

    void expose_emass()
    {
      namespace bp = boost::python;
      bp::class_< eMass >( "eMass", "Functor for computing effective masses" )
        .def_readwrite( "order", &eMass::order, "Order of interpolation." )
        .def_readwrite( "npoints", &eMass::npoints, "Number of points in interpolation." )
        .def_readwrite( "stepsize", &eMass::stepsize, "Distance between interpolation points." )
        .def_readwrite( "convergence", &eMass::tolerance, "Convergence criteria for cgs." )
        .def_readwrite( "verbose", &eMass::verbose, "Verbosity of cgs." )
        .def_readwrite( "itermax", &eMass::itermax, "Maximum number of iterations for cgs." )
        .add_property
         ( 
           "kpoint",
           bp::make_function(&eMass::get_kpoint, bp::return_value_policy<bp::return_by_value>()),
          &eMass::set_kpoint,
           "Kpoint at which to compute emass." 
         )
        .add_property
         ( 
           "direction",
           bp::make_function(&eMass::get_direction, bp::return_value_policy<bp::return_by_value>()),
          &eMass::set_direction,
           "Direction of emass." 
         )
        .def_readwrite( "nbstates", &eMass::nbstates, "Number of emasses to compute." )
        .def
        ( 
          "__call__", &get_masses2,
          (
            bp::arg("self"),
            bp::arg("escan"),
            bp::arg("ocell"),
            bp::arg("structure"),
            bp::arg("ref")
          ) 
        ) 
        .def
        ( 
          "__call__", &get_masses,
          (
            bp::arg("self"),
            bp::arg("escan"),
            bp::arg("ocell"),
            bp::arg("structure"),
            bp::arg("ref")
          ),
          "Computes and returns tuples ( eigenvalue at kpoint, effective mass at kpoint ).\n\n"
          "@param escan: the escan functional for computing eigenvalues.\n"
          "@type escan: L{Escan}\n"
          "@param ocell: original (unrelaxed) structure.\n"
          "@type ocell: 3x3 numpy array\n"
          "@param structure: relaxed structure.\n"
          "@type structure: L{lada.crystal.Structure}\n"
          "@param ref: reference energy at which to compute bands.\n"
        );
    }

  } // namespace python
} // namespace LaDa
