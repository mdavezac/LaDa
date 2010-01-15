//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/errors.hpp>
#include <sstream>

#include <opt/types.h>
#include <opt/debug.h>
#include <python/misc.hpp>
#include <python/xml.hpp>
#include <python/std_vector.hpp>

#include "../structure.h"
#include "../atom.h"
#include "../lattice.h"



#include "atom.hpp"


namespace LaDa
{
  namespace Python
  {
    template< class T_TYPE > Crystal::Atom_Type<T_TYPE>* default_constructor();
    template< class T_TYPE > 
      Crystal::Atom_Type<T_TYPE>* copy_constructor( const Crystal::Atom_Type<T_TYPE>& _ob );
    template< class T_TYPE >
      Crystal::Atom_Type<T_TYPE>* object_constructor( const boost::python::tuple& _ob );
    types::t_real toReal(std::string _str );
    std::string toType( types::t_real _r );

    template< class T_TYPE >
      void expose_typed_atom( const std::string &_name,
                              const std::string &_ds,
                              const std::string &_typeds )
      {
        namespace bp = boost::python;
        typedef Crystal::Atom_Type< T_TYPE > t_Atom;
        bp::class_< t_Atom >( _name.c_str(), _ds.c_str() )
          .def( "__init__", bp::make_constructor( default_constructor< T_TYPE > ) )
          .def( "__init__", bp::make_constructor( copy_constructor< T_TYPE > ) )
          .def( "__init__", bp::make_constructor( object_constructor< T_TYPE > ) )
          .def_readwrite( "pos",    &t_Atom::pos,
                          "a LaDa.rVector3d object containing the"
                          " atomic position in cartesian units." )
          .def_readwrite( "site",   &t_Atom::site,
                          "index of the \"site\" as referenced by a LaDa.Lattice object." )
          .def_readwrite( "type",   &t_Atom::type, _typeds.c_str() )
          .def_readwrite( "freeze", &t_Atom::freeze )
          .def( "__str__",  &print<t_Atom> );

        expose_vector< t_Atom >( ( _name + "s" ).c_str(), ("A list of " + _name ).c_str() );
      }

    void expose_atom()
    {
      namespace bp = boost::python;
      typedef Crystal::Atom_Type<std::string> t_StrAtom;
      typedef Crystal::Atom_Type<types::t_real> t_Atom;
      typedef Crystal::Structure::t_Atom t_Atom;
      typedef Crystal::Structure::t_kAtom t_kAtom;
      typedef Crystal::Lattice::t_Site t_Site;

      bp::enum_<t_Atom::t_FreezeAtom>( "FreezeAtom", "Tags to freeze atomic coordinates." )
        .value(       "none", t_Atom::FREEZE_NONE )
        .value(          "x", t_Atom::FREEZE_X )
        .value(          "y", t_Atom::FREEZE_Y )
        .value(          "z", t_Atom::FREEZE_Z )
        .value(       "type", t_Atom::FREEZE_T )
        .value( "cartesians", t_Atom::FREEZE_CARTESIANS )
        .value(        "all", t_Atom::FREEZE_ALL )
        .export_values();

      expose_typed_atom< t_Atom :: t_Type >
      (
        "Atom", 
        "Atom for which the type is specified as a real number",
        "Atomic specie as a real number (usually -1.0 and 1.0)."
      );
      expose_typed_atom< t_StrAtom :: t_Type >
      (
        "StrAtom", 
        "Atom for which the type is specified as a string",
        "Atomic specie as a string."
      );
      expose_typed_atom< t_kAtom :: t_Type >
      (
        "kAtom",
        "Represents a reciprocal-space vector.",
        "Complex real-value representing the intensity of this k-vector."
      );
      expose_typed_atom< t_Site :: t_Type >
      (
        "Site",
        "A lattice-site listing all possible atomic-specie occupation",
        "A list of all possible atomic-species (as atomic simbols)"
      );
      expose_vector< t_Site::t_Type::value_type >( "StringVector", "A std::vector of strings." );
      

      bp::def( "toAtomType", &toReal, "Converts from an atomic symbol to a real value." );
      bp::def( "fromAtomType", &toType, "Converts from a real value to an atomic symbol." );
    }

    void construct_type( Crystal::Atom_Type<types::t_real>& _atm,
                         const boost::python::tuple& _object )
    {
      const boost::python::object o( _object[3] );
      try{ _atm.type = boost::python::extract< types::t_real >( o ); }
      catch(...)
      {
        if( not Crystal::Structure::lattice )
        {
          PyErr_SetString(PyExc_RuntimeError,"Did you forget to initialize the Lattice?" );
          boost::python::throw_error_already_set();
          return;
        }
        Crystal::StrAtom stratom;
        stratom.pos = _atm.pos;
        stratom.type = boost::python::extract< std::string >( _object );
        Crystal::Structure::lattice->convert_StrAtom_to_Atom( stratom, _atm );
      }
    }
    void construct_type( Crystal::Atom_Type<types::t_complex>& _atm,
                         const boost::python::tuple& _object )
    {
      const boost::python::tuple o( _object[3] );
      _atm.type = types::t_complex( boost::python::extract< types::t_real >( o[0] ),
                                    boost::python::extract< types::t_real >( o[1] ) );
    }
    void construct_type( Crystal::Atom_Type<std::string>& _atm,
                         const boost::python::object& _object )
    {
      const boost::python::str o( _object[3] );
      _atm.type = boost::python::extract< std::string >( o ); 
    }
    void construct_type( Crystal::Atom_Type< std::vector<std::string> >& _atm,
                         const boost::python::tuple& _object )
    {
      const boost::python::list o( _object[3] );
      for( size_t i(0), n( boost::python::len( o ) ); i < n; ++n )
        _atm.type.push_back( boost::python::extract<std::string>( o[i]) );
    }

    template< class T_TYPE > 
      Crystal::Atom_Type<T_TYPE>* default_constructor()
        { return new Crystal::Atom_Type<T_TYPE>; }
    template< class T_TYPE >
      Crystal::Atom_Type<T_TYPE>* copy_constructor( const Crystal::Atom_Type<T_TYPE>& _ob )
        { return new Crystal::Atom_Type<T_TYPE>( _ob ); }
    template< class T_TYPE >
      Crystal::Atom_Type<T_TYPE>* object_constructor( const boost::python::tuple& _ob )
      {
        using namespace boost::python;
        typedef Crystal::Atom_Type<T_TYPE> t_Atom;
        t_Atom *result = NULL;
        try
        { 
          result = new t_Atom;
          types::t_unsigned length = len(_ob);
          if( length < 3 ) return result;
        
          result->pos.x[0] = extract< types::t_real >( _ob[0] );
          result->pos.x[1] = extract< types::t_real >( _ob[1] );
          result->pos.x[2] = extract< types::t_real >( _ob[2] );
          if( length == 3 ) return result;
          construct_type( *result, _ob );
          if( length == 4 ) return result;
          result->site = extract< types::t_int >( _ob[4] );
          return result;
        }
        catch( std::exception &_e )
        {
          if( result ) delete result;
          std::ostringstream sstr;
          sstr << "Object cannot be converted to an atom: \n"
               << _e.what() << "\n";
          PyErr_SetString(PyExc_RuntimeError, sstr.str().c_str());
          boost::python::throw_error_already_set();
        }
        catch( ... )
        {
          if( result ) delete result;
          PyErr_SetString(PyExc_RuntimeError, "Could not convert object to Atom." );
          boost::python::throw_error_already_set();
        }
        return NULL;
      }

    types::t_real toReal(std::string _str )
    { 
      if( not Crystal::Structure::lattice )
      {
        PyErr_SetString(PyExc_RuntimeError, "Could not convert atom type.\n" );
        boost::python::throw_error_already_set();
        return -1e0;
      }
      if( _str.compare( Crystal::Structure::lattice->sites[0].type[0] ) == 0 )
        return types::t_real(-1);
      else if( _str.compare( Crystal::Structure::lattice->sites[0].type[1] ) == 0 )
        return types::t_real(1);
      else
      {
        PyErr_SetString(PyExc_RuntimeError, "Requested Atomic type is not within lattice" );
        boost::python::throw_error_already_set();
        return -1e0;
      }
    }
    std::string toType( types::t_real _r )
    { 
      if( not Crystal::Structure::lattice )
      {
        PyErr_SetString(PyExc_RuntimeError, "Could not convert atom type.\n" );
        boost::python::throw_error_already_set();
        return "-1";
      }
      if( Crystal::Structure::lattice->sites.size() != 1)
      {
        PyErr_SetString(PyExc_RuntimeError, "Lattice cannot have more than one atomic site.\n" ); 
        boost::python::throw_error_already_set();
        return "-1";
      }
      if( Crystal::Structure::lattice->sites[0].type.size() != 2)
      {
        PyErr_SetString(PyExc_RuntimeError, "Lattice must have two atomic types.\n" ); 
        boost::python::throw_error_already_set();
        return "-1";
      }
      if(     math::neq( _r, types::t_real(1) ) 
          and math::neq( _r, types::t_real(-1) ) ) return "";
      return math::neq( _r, types::t_real(1) ) ? 
                Crystal::Structure::lattice->sites[0].type[0]:
                Crystal::Structure::lattice->sites[0].type[1] ;
    }

  }
} // namespace LaDa
