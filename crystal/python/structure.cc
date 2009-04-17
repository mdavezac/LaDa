//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <sstream>
#include <complex>

#include <boost/python.hpp>
#ifdef _MPI
# include <boost/mpi/python.hpp>
#endif
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <python/misc.hpp>
#include <python/xml.hpp>
#include <atat/serialize.h>

#include "../structure.h"
#include "../read_poscar.h"
#include "../lattice.h"
#include "../fractional_cartesian.h"

#include "structure.hpp"


namespace LaDa
{
  namespace Python
  {
    namespace XML
    {
      template<> std::string nodename<Crystal::Structure>()  { return "Structure"; }
      template<> std::string nodename< Crystal::TStructure<std::string> >() 
        { return nodename<Crystal::Structure>(); }
    }

    template< class BULL >
    Crystal::Lattice& return_crystal_lattice( BULL & )
    {
      __DOASSERT( not Crystal::Structure::lattice,
                  "Lattice pointer has not been set.\n" )
      return *Crystal::Structure::lattice; 
    }

    template<class T_TYPE> 
      void read_poscar( Crystal::TStructure<T_TYPE> &_struct, 
                        const boost::python::str& _path,
                        const boost::python::list& _types )
      {
        const std::string str = boost::python::extract<std::string>( _path );
        const boost::filesystem::path path( str );
        namespace bp = boost::python;
        const size_t nb( bp::len( _types ) );
        __DOASSERT( nb == 0, "No types given on input to read_poscar.\n" )
        std::vector< T_TYPE > types(nb);
        for( size_t i(0); i < nb; ++i )
          types[i] = bp::extract<T_TYPE>( _types[i] );
        Crystal :: read_poscar< T_TYPE >( _struct, path, types ); 
      }

    template< class T_STRUCTURE >
      struct pickle_structure : boost::python::pickle_suite
      {
        static boost::python::tuple getinitargs( T_STRUCTURE const& _w)  
        {
          return boost::python::tuple();
        }
        static boost::python::tuple getstate(const T_STRUCTURE& _in)
        {
          std::ostringstream ss;
          boost::archive::text_oarchive oa( ss );
          oa << _in;

          return boost::python::make_tuple( ss.str() );
        }
        static void setstate( T_STRUCTURE& _out, boost::python::tuple state)
        {
          namespace bp = boost::python;
          if( bp::len( state ) != 1 )
          {
            PyErr_SetObject(PyExc_ValueError,
                            ("expected 1-item tuple in call to __setstate__; got %s"
                             % state).ptr()
                );
            bp::throw_error_already_set();
          }
          const std::string str = bp::extract< std::string >( state[0] );
          std::istringstream ss( str.c_str() );
          boost::archive::text_iarchive ia( ss );
          ia >> _out;
        }
      };

    void expose_structure()
    {
      namespace bp = boost::python;
      bp::class_< Crystal::Structure >( "Structure", "Defines a structure.\n"
                                        "Generally, it is a super-cell of a LaDa.Lattice object." )
        .def( bp::init< Crystal::Structure& >() )
        .def_readwrite( "cell",    &Crystal::Structure::cell,
                        "The cell in cartesian coordinates (in units of LaDa.Structure.scale)." )
        .def_readwrite( "atoms",   &Crystal::Structure::atoms,
                        "The list of atoms of type LaDa.details_Atom. "
                        "Coordinates are in units of LaDa.Structure.Scale" )
        .def_readwrite( "k_vecs",  &Crystal::Structure::k_vecs,
                        "The list of reciprocal-space vectors."
                        " It is constructure with respected to a LaDa.Lattice object.\n"  ) 
        .def_readwrite( "energy",  &Crystal::Structure::energy, "Holds a real value." )
        .def_readwrite( "weight",  &Crystal::Structure::weight, "Optional weight for fitting purposes." )
        .def_readwrite( "scale",   &Crystal::Structure::scale,
                        "A scaling factor for atomic-positions and cell-vectors." )
        .def_readwrite( "name", &Crystal::Structure::name, "Holds a string.\n" )
        .def( "__str__",  &print<Crystal::Structure> ) 
        .def( "fromXML",  &XML::from<Crystal::Structure>, bp::arg("file"),
              "Loads a structure from an XML file." )
        .def( "toXML",  &XML::to<Crystal::Structure>, bp::arg("file"),
              "Adds a tag to an XML file describing this structure."  )
        .def( "lattice", &return_crystal_lattice< Crystal::Structure >,
              bp::return_value_policy<bp::reference_existing_object>(),
              "References the lattice within which this structure is defined."
              " Read, but do not write to this object." )
        .def( "concentration", &Crystal::Structure::get_concentration, "Returns concentration." )
        .def_pickle( pickle_structure< Crystal::Structure >() );
      bp::class_< Crystal::TStructure<std::string> >( "sStructure" )
        .def( bp::init< Crystal::TStructure<std::string>& >() )
        .def_readwrite( "cell",    &Crystal::TStructure<std::string>::cell,
                        "The cell in cartesian coordinates (in units of LaDa.Structure.scale)." )
        .def_readwrite( "atoms",   &Crystal::TStructure<std::string>::atoms,
                        "The list of atoms of type LaDa.details_Atom. "
                        "Coordinates are in units of LaDa.Structure.Scale" )
        .def_readwrite( "energy",  &Crystal::TStructure<std::string>::energy, "Holds a real value." )
        .def_readwrite( "weight",  &Crystal::TStructure<std::string>::weight,
                        "Optional weight for fitting purposes." )
        .def_readwrite( "scale",   &Crystal::TStructure<std::string>::scale,
                        "A scaling factor for atomic-positions and cell-vectors." )
        .def_readwrite( "name", &Crystal::TStructure<std::string>::name, "Holds a string." )
        .def( "__str__",  &print< Crystal::TStructure<std::string> > ) 
        .def( "fromXML",  &XML::from< Crystal::TStructure<std::string> >, bp::arg("file"),
              "Loads a structure from an XML file." )
        .def( "toXML",  &XML::to< Crystal::TStructure<std::string> >, bp::arg("file"),
              "Adds a tag to an XML file describing this structure."  )
        .def( "lattice", &return_crystal_lattice< Crystal::TStructure<std::string> >,
              bp::return_value_policy<bp::reference_existing_object>(),
              "References the lattice within which this structure is defined."
              "Read, but do not write to this object." )
        .def_pickle( pickle_structure< Crystal::TStructure<std::string> >() );
      bp::def("read_poscar", &read_poscar<std::string>,
              (
                bp::arg("structure"),
                bp::arg("filename"),
                bp::arg("species")
              ),
              "Reads a vasp POSCAR and fills in structure object.\n"
              "Third argument is a list of types, since types are implicit in POSCAR files." );
      bp::def("to_cartesian", &Crystal::to_cartesian<std::string>,
              "Transforms a structure from cartesian to fractional coordinates.\n" );
      bp::def("to_fractional", &Crystal::to_fractional<std::string>,
              "Transforms a structure from fractional to cartesian coordinates.\n" );
    }

  }
} // namespace LaDa
