//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <sstream>
#include <complex>

#include <boost/shared_ptr.hpp>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/str.hpp>
#include <boost/python/other.hpp>
#include <boost/python/self.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/reference_existing_object.hpp>
#include <boost/python/register_ptr_to_python.hpp>
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
#include <math/serialize.h>

#include <physics/physics.h>

#include "../structure.h"
#include "../fill_structure.h"
#include "../lattice.h"
#include "../fractional_cartesian.h"

#include "structure.hpp"


namespace LaDa
{
  namespace Python
  {
    namespace bp = boost::python;
    namespace XML
    {
      template<> std::string nodename<Crystal::Structure>()  { return "Structure"; }
      template<> std::string nodename< Crystal::TStructure<std::string> >() 
        { return nodename<Crystal::Structure>(); }
    }

    template< class BULL >
    Crystal::Lattice const* return_crystal_lattice( BULL & )
    {
      if( not Crystal::Structure::lattice )
      {
        PyErr_SetString(PyExc_RuntimeError,
                        "Crystal::Structure::lattice has not been set.\n");
        bp::throw_error_already_set();
        return NULL;
      }
      return Crystal::Structure::lattice; 
    }

    template< class T_STRUCTURE >
      struct pickle_structure : bp::pickle_suite
      {
        static bp::tuple getinitargs( T_STRUCTURE const& _w)  
        {
          return bp::tuple();
        }
        static bp::tuple getstate(const T_STRUCTURE& _in)
        {
          std::ostringstream ss;
          boost::archive::text_oarchive oa( ss );
          oa << _in;

          return bp::make_tuple( ss.str() );
        }
        static void setstate( T_STRUCTURE& _out, bp::tuple state)
        {
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

    bp::str xcrysden( Crystal::Structure const & _str )
    {
      std::ostringstream sstr;
      _str.print_xcrysden( sstr ); 
      return bp::str( sstr.str().c_str() );
    }
    bp::str xcrysden_str( Crystal::TStructure<std::string> const & _str )
    {
      if( not _str.lattice ) return bp::str("");
      std::ostringstream sstr;
      sstr << "CRYSTAL\nPRIMVEC\n" << ( (~_str.cell) * _str.scale ) << "\nPRIMCOORD\n" 
                << _str.atoms.size() << " 1 \n";  
      Crystal::TStructure<std::string>::t_Atoms :: const_iterator i_atom = _str.atoms.begin();
      Crystal::TStructure<std::string>::t_Atoms :: const_iterator i_atom_end = _str.atoms.end();
      for(; i_atom != i_atom_end; ++i_atom )
      {
        sstr << " " << Physics::Atomic::Z(i_atom->type) 
             << " " << ( i_atom->pos * _str.scale ) << "\n";
      }
      return bp::str( sstr.str().c_str() );
    }

    template< class T_STRUCTURE >
      T_STRUCTURE* empty()
      {
        T_STRUCTURE* result = new T_STRUCTURE(); 
        if( result->lattice ) result->scale = result->lattice->scale;
        return result;
      }
    
    template< class T_STRUCTURE >
      T_STRUCTURE* copy( T_STRUCTURE const &_o )
      {
        T_STRUCTURE* result = new T_STRUCTURE( _o ); 
        if( result->lattice ) result->scale = result->lattice->scale;
        return result;
      }

    Crystal::TStructure<std::string>* real_to_string( Crystal::Structure const &_s )
    {
      Crystal::TStructure<std::string>* result = new Crystal::TStructure<std::string>(); 
      Crystal::convert_real_to_string_structure( _s, *result );     
      return result;
    }
    Crystal::Structure* string_to_real( Crystal::TStructure<std::string> const &_s )
    {
      Crystal::Structure* result = new Crystal::Structure(); 
      Crystal::convert_string_to_real_structure( _s, *result );     
      return result;
    }

    void expose_structure()
    {
      typedef Crystal::Structure::t_FreezeCell t_FreezeCell;
      bp::enum_<t_FreezeCell>( "FreezeCell", "Tags to freeze cell coordinates." )
        .value( "none", Crystal::Structure::FREEZE_NONE )
        .value(   "xx", Crystal::Structure::FREEZE_XX )
        .value(   "xy", Crystal::Structure::FREEZE_XY )
        .value(   "xz", Crystal::Structure::FREEZE_XZ )
        .value(   "yy", Crystal::Structure::FREEZE_YY )
        .value(   "yz", Crystal::Structure::FREEZE_YZ )
        .value(   "zz", Crystal::Structure::FREEZE_ZZ )
        .value(  "all", Crystal::Structure::FREEZE_ALL )
        .export_values();

      bp::class_< Crystal::Structure >( "Structure", "Defines a structure.\n"
                                        "Generally, it is a super-cell of a LaDa.Lattice object." )
        .def( "__init__", bp::make_constructor( copy<Crystal::Structure> ) )
        .def( "__init__", bp::make_constructor( empty<Crystal::Structure> ) )
        .def( "__init__", bp::make_constructor( string_to_real ) )
        .def_readwrite( "cell",    &Crystal::Structure::cell,
                        "The cell in cartesian coordinates (in units of LaDa.Structure.scale)." )
        .def_readwrite( "freeze", &Crystal::Structure::freeze,
                        "Tags to freeze coordinates when relaxing structure.\n" )
        .def_readwrite( "atoms",   &Crystal::Structure::atoms,
                        "The list of atoms of type LaDa.details_Atom. "
                        "Coordinates are in units of LaDa.Structure.Scale" )
        .def_readwrite( "k_vecs",  &Crystal::Structure::k_vecs,
                        "The list of reciprocal-space vectors."
                        " It is constructure with respected to a LaDa.Lattice object.\n"  ) 
        .def_readwrite( "energy",  &Crystal::Structure::energy, "Holds a real value." )
        .def_readwrite( "weight",  &Crystal::Structure::weight,
                        "Optional weight for fitting purposes." )
        .def_readwrite( "scale",   &Crystal::Structure::scale,
                        "A scaling factor for atomic-positions and cell-vectors." )
        .def_readwrite( "name", &Crystal::Structure::name, "Holds a string.\n" )
        .def( "__str__",  &print<Crystal::Structure> ) 
        .def( "fromXML",  &XML::from<Crystal::Structure>, bp::arg("file"),
              "Loads a structure from an XML file." )
        .def( "toXML",  &XML::to<Crystal::Structure>, bp::arg("file"),
              "Adds a tag to an XML file describing this structure."  )
        .add_property
        ( 
          "lattice",
          bp::make_function
          (
            &return_crystal_lattice< Crystal::Structure >,
            bp::return_value_policy<bp::reference_existing_object>()
          ),
          "References the lattice within which this structure is defined."
          " Read, but do not write to this object." 
        )
        .def( "concentration", &Crystal::Structure::get_concentration, "Returns concentration." )
        .def( bp::self == bp::other<Crystal::Structure>() )
        .def( "xcrysden", &xcrysden, "Outputs in XCrysden format." )
        .def_pickle( pickle_structure< Crystal::Structure >() );
      bp::class_< Crystal::TStructure<std::string> >( "sStructure" )
        .def( "__init__", bp::make_constructor( copy< Crystal::TStructure<std::string> > ) )
        .def( "__init__", bp::make_constructor( empty< Crystal::TStructure<std::string> > ) )
        .def( "__init__", bp::make_constructor( real_to_string ) )
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
        .def( "xcrysden", &xcrysden_str, "Outputs in XCrysden format." )
        .add_property
        ( 
          "lattice",
          bp::make_function
          (
            &return_crystal_lattice< Crystal::TStructure<std::string> >,
            bp::return_value_policy<bp::reference_existing_object>()
          ),
          "References the lattice within which this structure is defined."
          " Read, but do not write to this object." 
        )
        .def_pickle( pickle_structure< Crystal::TStructure<std::string> >() );
      bp::def("to_cartesian", &Crystal::to_cartesian<std::string>,
              "Transforms a structure from cartesian to fractional coordinates.\n" );
      bp::def("to_fractional", &Crystal::to_fractional<std::string>,
              "Transforms a structure from fractional to cartesian coordinates.\n" );

      bp::register_ptr_to_python< boost::shared_ptr<Crystal::Structure> >();
      bp::register_ptr_to_python< boost::shared_ptr< Crystal::TStructure<std::string> > >();

      bool (*for_real)( Crystal::Structure& ) = &Crystal::fill_structure;
      bool (*for_string)( Crystal::TStructure<std::string>& ) = &Crystal::fill_structure;
      bp::def("fill_structure", for_real, "Fills a structure when atomic positions are unknown." );
      bp::def("fill_structure", for_string, "Fills a structure when atomic positions are unknown." );
    }

  }
} // namespace LaDa
