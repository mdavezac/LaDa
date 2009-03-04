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
    boost::python::list get_masses
    (
      const Pescan::eMass& _emass,
      const Pescan::Interface& _interface,
      const atat::rMatrix3d &_ocell, 
      const Crystal::Structure &_structure,
      const atat::rVector3d &_at,
      const atat::rVector3d &_direction,
      const size_t &_nbstates,
      const types::t_real &_eref
    )
    {
      typedef std::vector< std::pair< types::t_real, types::t_real > > t_Out;
      t_Out out;
      _emass
      (
       _interface,
       _ocell, 
       _structure,
       _at,
       _direction,
       _nbstates,
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
      typedef LaDa::Pescan::eMass t_eMass;
      bp::class_< t_eMass >( "eMass", "Functor for computing effective masses" )
        .def_readwrite( "cgs", &t_eMass::cgs, "Conjugate gradient for least-square-fit." ) 
        .def_readwrite( "order", &t_eMass::order, "Order of interpolation." )
        .def_readwrite( "npoints", &t_eMass::npoints, "Number of points in interpolation." )
        .def_readwrite( "stepsize", &t_eMass::stepsize, "Distance between interpolation points." )
        .def
        ( 
          "__call__", &get_masses,
          (
            bp::arg("self"),
            bp::arg("escan"),
            bp::arg("ocell"),
            bp::arg("structure"),
            bp::arg("kpoint") = atat::rVector3d(0,0,0),
            bp::arg("direction"),
            bp::arg("nbstates") = 2,
            bp::arg("ref")
          ),
          "Computes and returns tuples ( eigenvalue at kpoint, effective mass at kpoint ).\n"
          "escan = the escan functional for computing eigenvalues.\n"
          "ocell = original (unrelaxed) structure.\n"
          "structure = relaxed structure.\n"
          "kpoint = reciprocal space vector for which to compute effective masses.\n"
          "direction = direction for which compute effective masses.\n"
          "nbstates = number of bands for which compute effective masses. Should be even.\n"
          "ref = reference energy at which to compute bands.\n"
        );
    }

  } // namespace Python
} // namespace LaDa
