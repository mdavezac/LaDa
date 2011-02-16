#include "LaDaConfig.h"

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
        "Index of the site which matches the given position.\n\n"
        ":Parameters:\n"\
        "  pos\n    A position to match modulo the supercell.\n"
        "  invcell\n    the inverse of the supercell.\n"
        "  sites : `Lattice.sites`, `Structure.atoms`, `rStructure.atoms`\n"
        "    A list of sites.\n"
        // other implementation,
        "  lattice : `Lattice`\n"
        "    match pos to the sites in lattice, modulo the lattice vectors.\n"
        // other implementation,
        "  structure : `Structure`\n"
        "    match pos to the atoms in structure, modulo the cell-vectors.\n\n"
        "This function can take either three (pos, invcell, sites) or two "
        "arguments (pos, structure), or (pos, lattice).\n"
      );
    }

  }
} // namespace LaDa
