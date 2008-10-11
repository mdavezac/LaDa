//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>
#include <opt/types.h>
#include <crystal/structure.h>
#include <crystal/lattice.h>

#include "xml.hpp"
#include "misc.hpp"
#include "lattice.hpp"

namespace PythonLaDa
{
  namespace XML
  {
    template<> std::string nodename<Crystal::Lattice>()    { return "Lattice"; }
    template<> void do_specialcode< Crystal::Lattice >( Crystal::Lattice &_type )
    {
      Crystal::Structure::lattice = &_type; 
      _type.find_space_group();
    }

  }

  boost::python::list symmetry_ops( const Crystal::Lattice& _lat )
  {
    types::t_int N( _lat.space_group.point_op.getSize() );
    boost::python::list result;
    for( types::t_int i = 0; i < N; ++i )
    {
      boost::python::list inner;
      inner.append( _lat.space_group.point_op(i) );
      inner.append( _lat.space_group.trans(i) );
      result.append( inner );
    }
    return result;
  }

// Crystal::Lattice& return_crystal_lattice()
// {
//   __DOASSERT( not Crystal::Structure::lattice,
//               "Lattice pointer has not been set.\n" )
//   return *Crystal::Structure::lattice; 
// }

  void expose_lattice()
  {
    using namespace boost::python;
    class_< Crystal::Lattice >( "Lattice" )
      .def( init< Crystal::Lattice >() )
      .def_readwrite( "cell",  &Crystal::Lattice::cell )
      .def_readwrite( "sites", &Crystal::Lattice::sites )
      .def_readwrite( "scale", &Crystal::Lattice::scale )
      .def( "__str__",  &print<Crystal::Lattice> )
      .def( "syms",  &symmetry_ops )
      .def( "fromXML",  &XML::from<Crystal::Lattice> );
  // def( "StructureLattice", &return_crystal_lattice, 
  //      return_value_policy<reference_existing_object>() );
  }

}
