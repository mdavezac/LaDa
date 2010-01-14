//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/docstring_options.hpp>

#include <crystal/lattice.h>
#include <crystal/structure.h>

#include "../which_site.h"
#include "which_site.hpp"


namespace LaDa
{
  namespace Python
  {
    types::t_int which_site1( atat::rVector3d const &_pos, Crystal::Lattice const &_lat )
      { return Crystal::which_site( _pos, !_lat.cell, _lat.sites ); }
    types::t_int which_site2( atat::rVector3d const &_pos, Crystal::Structure const &_lat )
      { return Crystal::which_site( _pos, !_lat.cell, _lat.atoms ); }
    types::t_int which_site3( atat::rVector3d const &_pos,
                              Crystal::TStructure<std::string> const &_lat )
      { return Crystal::which_site( _pos, !_lat.cell, _lat.atoms ); }
    void expose_which_site()
    {
      boost::python::def
      ( 
        "which_site",
        &Crystal::which_site<Crystal::Structure::t_Atoms>,
        ( boost::python::arg("pos"), boost::python::arg("invcell"), boost::python::arg("sites") )
      );
      boost::python::def
      ( 
        "which_site",
        &Crystal::which_site<Crystal::TStructure<std::string>::t_Atoms>,
        ( boost::python::arg("pos"), boost::python::arg("invcell"), boost::python::arg("sites") )
      );
      boost::python::def
      ( 
        "which_site",
        &which_site1,
        ( boost::python::arg("pos"), boost::python::arg("lattice") )
      );
      boost::python::def
      ( 
        "which_site",
        &which_site2,
        ( boost::python::arg("pos"), boost::python::arg("structure") )
      );
      boost::python::def
      ( 
        "which_site",
        &which_site3,
        ( boost::python::arg("pos"), boost::python::arg("structure") )
      );
      boost::python::def
      ( 
        "which_site",
        &Crystal::which_site<Crystal::Lattice::t_Sites>,
        ( boost::python::arg("pos"), boost::python::arg("invcell"), boost::python::arg("sites") ),
        "Returns index of the site which matches the given position.\n\n"
        "This function can take either three (pos, invcell, sites) or two "
        "arguments (pos, structure), or (pos, lattice).\n"
        "@param pos: A position to match modulo the supercell.\n"
        "@type pos: L{atat.rVector3d}\n"
        "@param invcell: the inverse of the supercell.\n"
        "@type invcell: L{atat.rMatrix3d}\n"
        "@param sites: A list of sites.\n"
        "@type sites: L{Lattice.sites}, L{Structure.atoms}, L{sStructure.atoms}\n"
        // other implementation,
        "@param lattice: match pos to the sites in lattice, modulo the lattice vectors.\n"
        "@type lattice: L{Lattice}\n"
        // other implementation,
        "@param structure: match pos to the atoms in structure, modulo the cell-vectors.\n"
        "@type structure: L{Structure} or L{sStructure}\n"
      );
    }

  }
} // namespace LaDa
