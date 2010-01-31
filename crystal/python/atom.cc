//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
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
    template< class T_TYPE >
      Crystal::Atom_Type<T_TYPE>* object_constructor( math::rVector3d const &_vec,
                                                      boost::python::object const & _ob) ;
    template< class T_TYPE >
      Crystal::Atom_Type<T_TYPE>* full_constructor( math::rVector3d const &_vec,
                                                    boost::python::object const & _ob,
                                                    types::t_int _site );
    types::t_real toReal(std::string _str );
    std::string toType( types::t_real _r );

    template< class T_TYPE >
      void expose_typed_atom( const std::string &_name,
                              const std::string &_ds,
                              const std::string &_typeds )
      {
        namespace bp = boost::python;
        typedef Crystal::Atom_Type< T_TYPE > t_Atom;
        bp::class_< t_Atom >
        ( 
          _name.c_str(), 
          ( 
              _ds 
            +  "\nThis object can be constructed from:\n"
               "  - with no argument\n"
               "  - another L{" + _name + "} object (deepcopy)\n"
               "  - a numpy 1d-vector of length 3, in which self.L{pos<" 
               + _name + ".pos>} is the only variable to be set.\n"
               "  - same as above + a type value.\n"
               "  - same as above + a site index.\n"
          ).c_str()
        ).def(bp::init<t_Atom const &>())
         .def("__init__", bp::make_constructor( &object_constructor< T_TYPE > ) )
         .def("__init__", bp::make_constructor( &full_constructor< T_TYPE > ) )
         .def(bp::init<math::rVector3d const &, T_TYPE>())
         .add_property
         (
           "pos",
           make_getter(&t_Atom::pos, bp::return_value_policy<bp::return_by_value>()),
           make_setter(&t_Atom::pos, bp::return_value_policy<bp::return_by_value>()),
           "A 1-dimensional numpy array of length 3 containing atomic position in cartesian units."
         )
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
                         const boost::python::object& _object )
    {
      try{ _atm.type = boost::python::extract< types::t_real >( _object ); }
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
                         const boost::python::object& _object )
    {
      if( bp::len(_object) == 2 )
      {
        PyErr_SetString( PyExc_TypeError, "Object is not a complex value.\n");
        bp::throw_error_already_set();
        return; 
      }
      _atm.type = types::t_complex( boost::python::extract< types::t_real >( _object[0] ),
                                    boost::python::extract< types::t_real >( _object[1] ) );
    }
    void construct_type( Crystal::Atom_Type<std::string>& _atm,
                         const boost::python::object& _object )
      { _atm.type = boost::python::extract< std::string >( _object ); }
    void construct_type( Crystal::Atom_Type< std::vector<std::string> >& _atm,
                         const boost::python::object& _object )
    {
      namespace bp = boost::python;
      PyObject * const obj_ptr( _object.ptr() );
      if( PyString_Check(obj_ptr) )
      {
        _atm.type.push_back( bp::extract<std::string>(_object) );
        return;
      }
      else if( PySequence_Check(obj_ptr) )
      {
        size_t const N = PySequence_Length(obj_ptr);
        for( size_t i(0); i < N; ++i )
        {
          PyObject * const item_ptr( PySequence_GetItem(obj_ptr, i) );
          if( not PyString_Check(item_ptr) )
          {
            PyErr_SetString(PyExc_TypeError, "Object in input sequence is not a string.\n");
            bp::throw_error_already_set();
            return;
          }
          _atm.type.push_back( PyString_AsString(item_ptr) );
        }
      }
      else
      {
        PyErr_SetString
        (
          PyExc_TypeError, 
          "Input object is neither a string nor a sequence of strings.\n"
        );
        bp::throw_error_already_set();
        return;
      }
    }

    template< class T_TYPE >
      Crystal::Atom_Type<T_TYPE>* object_constructor( math::rVector3d const &_vec,
                                                      boost::python::object const& _ob )
      {
        using namespace boost::python;
        typedef Crystal::Atom_Type<T_TYPE> t_Atom;
        t_Atom *result = NULL;
        try
        { 
          result = new t_Atom;
          result->pos = _vec;
          construct_type( *result, _ob );
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
    template< class T_TYPE >
      Crystal::Atom_Type<T_TYPE>* full_constructor( math::rVector3d const &_vec,
                                                    boost::python::object const& _ob,
                                                    types::t_int _site )
      {
        using namespace boost::python;
        typedef Crystal::Atom_Type<T_TYPE> t_Atom;
        t_Atom *result = NULL;
        try
        { 
          result = new t_Atom;
          result->pos = _vec;
          construct_type( *result, _ob );
          result->site = _site; 
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
