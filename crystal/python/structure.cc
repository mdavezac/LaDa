#include "LaDaConfig.h"

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
#include <boost/python/return_by_value.hpp>
#include <boost/python/reference_existing_object.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#ifdef LADA_MPI
# include <boost/mpi/python.hpp>
#endif
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/filesystem/operations.hpp>

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
             << " " << ( i_atom->pos[0] * _str.scale ) << " "
             << " " << ( i_atom->pos[1] * _str.scale ) << " "
             << " " << ( i_atom->pos[2] * _str.scale ) << "\n";
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

    template<class T_STRUCTURE> 
      T_STRUCTURE* fromXML(std::string const &_path)
      {
        namespace bfs = boost::filesystem;
        if( not bfs::exists(_path) )
        {
          PyErr_SetString(PyExc_IOError, (_path + " does not exist.\n").c_str());
          bp::throw_error_already_set();
          return NULL;
        }
        T_STRUCTURE* result = new T_STRUCTURE;
        try
        { 
          if(not result->Load(_path))
          {
            PyErr_SetString(PyExc_IOError, ("Could not load structure from " + _path).c_str());
            delete result;
            result = NULL;
            bp::throw_error_already_set();
          }
        }
        catch(std::exception &_e)
        {
          PyErr_SetString(PyExc_IOError, ("Could not load structure from " + _path).c_str());
          delete result;
          result = NULL;
          bp::throw_error_already_set();
        }
        return result;
      }

    template<class T_STRUCTURE>
      boost::shared_ptr<T_STRUCTURE> fill_structure(T_STRUCTURE const &_str)
      {
        boost::shared_ptr<T_STRUCTURE> result( new T_STRUCTURE(_str) );
        if( not result->lattice ) 
        {
          PyErr_SetString(PyExc_RuntimeError, "lada.crystal.lattice not set.\n");
          bp::throw_error_already_set();
          boost::shared_ptr<T_STRUCTURE> b;
          result.swap(b);
        }
        else if( not Crystal::fill_structure(*result) )
        {
          PyErr_SetString(PyExc_RuntimeError, "Could not create fill in structure.\n");
          bp::throw_error_already_set();
          boost::shared_ptr<T_STRUCTURE> b;
          result.swap(b);
        }
        return result;
      }
    boost::shared_ptr< Crystal::TStructure<std::string> >
      fill_cell(math::rMatrix3d const& _cell)
      {
        typedef Crystal::TStructure<std::string> t_structure;
        t_structure str;
        str.cell = _cell;
        boost::shared_ptr< t_structure > result(fill_structure<t_structure>(str));
        if( not result ) return result;
        result->scale = result->lattice->scale;
        return result;
      }
      

    template< class T>
      math::rMatrix3d get_cell( T const &_str ) { return _str.cell; }
    template< class T>
      void set_cell( T &_str, math::rMatrix3d const &_cell) {_str.cell = _cell; }

    template<class T_STRUCTURE>
      bp::class_<T_STRUCTURE>   expose( std::string const &_name,
                                        std::string const &_desc, 
                                        std::string const &_type ) 
      {
        return bp::class_<T_STRUCTURE>( _name.c_str(), _desc.c_str() )
          .def( bp::init<T_STRUCTURE const&>() )
          .def( "__init__", bp::make_constructor( fromXML<Crystal::Structure> ) )
          .add_property
          (
            "cell",
            bp::make_function
            (
              &get_cell<T_STRUCTURE>, 
              bp::return_value_policy<bp::return_by_value>()
            ), &set_cell<T_STRUCTURE>, 
            ("A 3x3 numpy array representing the cell vector in cartesian units, "
             "in units of self.L{scale<lada.crystal." + _name + ".scale>}.").c_str()
          )
          .def_readwrite( "atoms",   &T_STRUCTURE::atoms,
                          (   "The list of atoms of type L{" + _type
                            + "}, in units of self.{scale<" + _name + ">}.").c_str() )
          .def_readwrite( "energy",  &T_STRUCTURE::energy, "Holds a real value." )
          .def_readwrite( "weight",  &T_STRUCTURE::weight,
                          "Optional weight for fitting purposes." )
          .def_readwrite( "scale",   &T_STRUCTURE::scale,
                          "A scaling factor for atomic-positions and cell-vectors." )
          .def_readwrite( "name", &T_STRUCTURE::name, "Holds a string." )
          .def( "__str__",  &print< T_STRUCTURE > ) 
          .def( "fromXML",  &XML::from< T_STRUCTURE >, bp::arg("file"),
                "Loads a structure from an XML file." )
          .def( "toXML",  &XML::to< T_STRUCTURE >, bp::arg("file"),
                "Adds a tag to an XML file describing this structure."  )
          .def( "xcrysden", &xcrysden_str, "Outputs in XCrysden format." )
          .def_readwrite( "freeze", &T_STRUCTURE::freeze,
                           "Tags to freeze coordinates when relaxing structure.\n\n" 
                           "See L{FreezeCell} for possible values." 
                        )
          .def(bp::self == bp::other<T_STRUCTURE>())
          .add_property
          ( 
            "lattice",
            bp::make_function
            (
              &return_crystal_lattice< T_STRUCTURE >,
              bp::return_value_policy<bp::reference_existing_object>()
            ),
            "References the lattice within which this structure is defined."
            " Read, but do not write to this object." 
          ) 
          .def_pickle( Python::pickle< T_STRUCTURE >() );
      }

    void expose_structure()
    {
      typedef Crystal::Structure::t_FreezeCell t_FreezeCell;
      bp::enum_<t_FreezeCell>( "FreezeCell", "Tags to freeze cell coordinates." )
        .value( "none", Crystal::Structure::FREEZE_NONE )
        .value(   "xx", Crystal::Structure::FREEZE_XX )
        .value(   "xy", Crystal::Structure::FREEZE_XY )
        .value(   "xz", Crystal::Structure::FREEZE_XZ )
        .value(   "yx", Crystal::Structure::FREEZE_YX )
        .value(   "yy", Crystal::Structure::FREEZE_YY )
        .value(   "yz", Crystal::Structure::FREEZE_YZ )
        .value(   "zx", Crystal::Structure::FREEZE_ZX )
        .value(   "zy", Crystal::Structure::FREEZE_ZY )
        .value(   "zz", Crystal::Structure::FREEZE_ZZ )
        .value(  "all", Crystal::Structure::FREEZE_ALL )
        .value(  "a0", Crystal::Structure::FREEZE_A0 )
        .value(  "a1", Crystal::Structure::FREEZE_A1 )
        .value(  "a2", Crystal::Structure::FREEZE_A2 )
        .export_values();

      expose< Crystal::Structure >
      (
        "rStructure", 
        "Defines a structure.\n\nGenerally, it is a super-cell of a L{Lattice} object.",
        "rAtom"
      ).def( "__init__", bp::make_constructor( string_to_real ) )
       .def_readwrite( "k_vecs",  &Crystal::Structure::k_vecs,
                       "The list of reciprocal-space vectors."
                       " It is constructure with respected to a LaDa.Lattice object.\n"  ) 
       .def( "concentration", &Crystal::Structure::get_concentration, "Returns concentration." )
       .def( bp::self == bp::other<Crystal::Structure>() );
      bp::register_ptr_to_python< boost::shared_ptr<Crystal::Structure> >();
      expose< Crystal::TStructure<std::string> >
      (
        "Structure", 
        "Defines a structure.\n\nGenerally, it is a super-cell of a L{Lattice} object.",
        "Atom"
      ).def( "__init__", bp::make_constructor( real_to_string ) );
      bp::register_ptr_to_python< boost::shared_ptr< Crystal::TStructure<std::string> > >();

      bp::def("to_cartesian", &Crystal::to_cartesian<std::string>,
              "Transforms a structure from cartesian to fractional coordinates.\n" );
      bp::def("to_fractional", &Crystal::to_fractional<std::string>,
              "Transforms a structure from fractional to cartesian coordinates.\n" );

      bp::def("_fill_structure_impl", &fill_structure<Crystal::Structure>);
      bp::def("_fill_structure_impl", &fill_cell);
      bp::def
      (
        "_fill_structure_impl", 
        &fill_structure< Crystal::TStructure<std::string> >,
        "Returns a structure from knowledge of cell and lattice.\n\n"
        "The argument can be of type L{Structure}, L{rStructure}, "
        "or a numpy 3x3 float64 array. In the second case, the return is "
        "also a L{rStructure}. In all other cases, the return is an L{Structure}.\n"
        "@raise RuntimeError: If the filled structure could not be created.\n" 
      );
    }

  }
} // namespace LaDa
