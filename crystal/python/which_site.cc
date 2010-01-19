//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <crystal/lattice.h>
#include <crystal/structure.h>

#include "../which_site.h"
#include "which_site.hpp"


namespace LaDa
{
  namespace Python
  {
    types::t_int which_site1( math::rVector3d const &_pos, Crystal::Lattice const &_lat )
      { return Crystal::which_site( _pos, !_lat.cell, _lat.sites ); }
    types::t_int which_site2( math::rVector3d const &_pos, Crystal::Structure const &_lat )
      { return Crystal::which_site( _pos, !_lat.cell, _lat.atoms ); }
    types::t_int which_site3( math::rVector3d const &_pos,
                              Crystal::TStructure<std::string> const &_lat )
      { return Crystal::which_site( _pos, !_lat.cell, _lat.atoms ); }
    void expose_which_site()
    {
      boost::python::def
      ( 
        "which_site",
        &Crystal::which_site<Crystal::Lattice::t_Sites>,
        ( boost::python::arg("pos"), boost::python::arg("inv_cell"), boost::python::arg("sites") ),
        "Returns index of site in sites matching the position pos modulo !inv_cell."
      );
      boost::python::def
      ( 
        "which_site",
        &Crystal::which_site<Crystal::Structure::t_Atoms>,
        ( boost::python::arg("pos"), boost::python::arg("inv_cell"), boost::python::arg("sites") ),
        "Returns index of site in sites matching the position pos modulo !inv_cell."
      );
      boost::python::def
      ( 
        "which_site",
        &Crystal::which_site<Crystal::TStructure<std::string>::t_Atoms>,
        ( boost::python::arg("pos"), boost::python::arg("inv_cell"), boost::python::arg("sites") ),
        "Returns index of site in sites matching the position pos modulo !inv_cell."
      );
      boost::python::def
      ( 
        "which_site",
        &which_site1,
        ( boost::python::arg("pos"), boost::python::arg("lattice") ),
        "Returns index of site in lattice matching the position pos."
      );
      boost::python::def
      ( 
        "which_site",
        &which_site2,
        ( boost::python::arg("pos"), boost::python::arg("structure") ),
        "Returns index of site in structure matching the position pos."
      );
      boost::python::def
      ( 
        "which_site",
        &which_site3,
        ( boost::python::arg("pos"), boost::python::arg("structure") ),
        "Returns index of site in structure matching the position pos."
      );
    }

  }
} // namespace LaDa
