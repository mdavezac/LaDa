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
  namespace Python 
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


    void expose_emass()
    {
      namespace bp = boost::python;
      typedef eMass t_eMass;
      bp::class_< t_eMass >( "eMass", "Functor for computing effective masses" )
        .def_readwrite( "order", &t_eMass::order, "Order of interpolation." )
        .def_readwrite( "npoints", &t_eMass::npoints, "Number of points in interpolation." )
        .def_readwrite( "stepsize", &t_eMass::stepsize, "Distance between interpolation points." )
        .def_readwrite( "convergence", &t_eMass::tolerance, "Convergence criteria for cgs." )
        .def_readwrite( "verbose", &t_eMass::verbose, "Verbosity of cgs." )
        .def_readwrite( "itermax", &t_eMass::itermax, "Maximum number of iterations for cgs." )
        .def_readwrite( "kpoint", &t_eMass::kpoint, "Kpoint at which to compute emass." )
        .def_readwrite( "direction", &t_eMass::direction, "Direction of emass." )
        .def_readwrite( "nbstates", &t_eMass::nbstates, "Number of emasses to compute." )
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

  } // namespace Python
} // namespace LaDa
