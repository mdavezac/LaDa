//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <sstream>

#include <opt/types.h>
#include <opt/debug.h>
#include <crystal/structure.h>
#include <crystal/atom.h>


#include "misc.hpp"
#include "xml.hpp"

#include "atom.hpp"


namespace LaDa
{
  namespace Python
  {
    Crystal::TStructure<std::string>::t_Atom* SAtomFromObject( boost::python::list& _ob );
    Crystal::Structure::t_Atom* AtomFromObject( boost::python::list& _ob );
    Crystal::Lattice::t_Site* SiteFromObject( boost::python::list& _ob );

    types::t_real toReal(std::string _str );
    std::string toType( types::t_real _r );

    void expose_atom()
    {
      using namespace boost::python;
      typedef Crystal::TStructure<std::string>::t_Atom t_StrAtom;
      typedef Crystal::Structure::t_Atom t_Atom;
      typedef Crystal::Structure::t_kAtom t_kAtom;
      typedef Crystal::Lattice::t_Site t_Site;
      class_< Crystal::Structure::t_Atoms >("VecStrings")
        .def(vector_indexing_suite< t_Site::t_Type >());

      class_< t_StrAtom >( "details_StrAtom" )
        .def( init< t_StrAtom >() )
        .def_readwrite( "pos",    &t_StrAtom::pos )
        .def_readwrite( "site",   &t_StrAtom::site )
        .def_readwrite( "type",   &t_StrAtom::type )
        .def_readwrite( "freeze", &t_StrAtom::freeze )
        .def( "__str__",  &print<t_StrAtom> ) ;
      class_< t_Atom >( "details_Atom" )
        .def( init< t_Atom >() )
        .def_readwrite( "pos",    &t_Atom::pos )
        .def_readwrite( "site",   &t_Atom::site )
        .def_readwrite( "type",   &t_Atom::type )
        .def_readwrite( "freeze", &t_Atom::freeze )
        .def( "__str__",  &print<t_Atom> ) ;
      class_< t_kAtom >( "details_kAtom" )
        .def( init< t_kAtom >() )
        .def_readwrite( "pos",    &t_kAtom::pos )
        .def_readwrite( "site",   &t_kAtom::site )
        .def_readwrite( "type",   &t_kAtom::type )
        .def_readwrite( "freeze", &t_kAtom::freeze )
        .def( "__str__",  &print<t_kAtom> ) ;
      class_< t_Site >( "details_Site" )
        .def( init< t_Site >() )
        .def_readwrite( "pos",    &t_Site::pos )
        .def_readwrite( "site",   &t_Site::site )
        .def_readwrite( "type",   &t_Site::type )
        .def_readwrite( "freeze", &t_Site::freeze )
        .def( "__str__",  &print<t_Site> ) ;


      def( "SAtom", &SAtomFromObject,
           return_value_policy<manage_new_object>() );
      def( "Atom", &AtomFromObject,
           return_value_policy<manage_new_object>() );
      def( "Site", &SiteFromObject,
           return_value_policy<manage_new_object>() );

      class_< Crystal::TStructure<std::string>::t_Atoms >("SAtoms")
        .def(vector_indexing_suite< Crystal::TStructure<std::string>::t_Atoms >());
      class_< Crystal::Structure::t_Atoms >("Atoms")
        .def(vector_indexing_suite< Crystal::Structure::t_Atoms >());
      class_< Crystal::Structure::t_kAtoms >("kAtoms")
        .def(vector_indexing_suite< Crystal::Structure::t_kAtoms >());
      class_< Crystal::Lattice::t_Sites >("Sites")
        .def(vector_indexing_suite< Crystal::Lattice::t_Sites >());

      def( "toAtomType", &toReal );
      def( "fromAtomType", &toType );
    }

    Crystal::TStructure<std::string>::t_Atom* SAtomFromObject( boost::python::list& _ob )
    {
      using namespace boost::python;
      typedef Crystal::TStructure<std::string>::t_Atom t_Atom;
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
        result->pos = result->pos;
        result->type = extract< std::string >( _ob[3] );
        if( length == 4 ) return result;
        result->site = extract< types::t_real >( _ob[4] );
        return result;
      }
      catch( std::exception &_e )
      {
        if( result ) delete result;
        std::ostringstream sstr;
        sstr << "Object cannot be converted to an atom: \n"
             << _e.what() << "\n";
        throw std::runtime_error( sstr.str() );
      }
      catch( ... )
      {
        if( result ) delete result;
        throw std::runtime_error( "Could not convert object to Atom." );
      }
      return NULL;
    }

    Crystal::Structure::t_Atom* AtomFromObject( boost::python::list& _ob )
    {
      using namespace boost::python;
      typedef Crystal::Structure::t_Atom t_Atom;
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
        try{ result->type = extract< types::t_real >(_ob[3] ); }
        catch(...)
        {
          if( not Crystal::Structure::lattice )
            throw std::runtime_error( "Did you forget to initialize the Lattice?" );
          Crystal::StrAtom stratom;
          stratom.pos = result->pos;
          stratom.type = extract< std::string >( _ob[3] );
          Crystal::Structure::lattice->convert_StrAtom_to_Atom( stratom, *result );
        }
        if( length == 4 ) return result;
        result->site = extract< types::t_real >( _ob[4] );
        return result;
      }
      catch( std::exception &_e )
      {
        if( result ) delete result;
        std::ostringstream sstr;
        sstr << "Object cannot be converted to an atom: \n"
             << _e.what() << "\n";
        throw std::runtime_error( sstr.str() );
      }
      catch( ... )
      {
        if( result ) delete result;
        throw std::runtime_error( "Could not convert object to Atom." );
      }
      return NULL;
    }
    Crystal::Lattice::t_Site* SiteFromObject( boost::python::list& _ob )
    {
      using namespace boost::python;
      typedef Crystal::Lattice::t_Site t_Site;
      types::t_unsigned length = len( _ob );
      if( length < 3 )
        throw std::runtime_error( "Object cannot be converted to an atom" );

      t_Site *result = new t_Site;
      try
      { 
        result->pos.x[0] = extract< types::t_real >( _ob[0] );
        result->pos.x[1] = extract< types::t_real >( _ob[1] );
        result->pos.x[2] = extract< types::t_real >( _ob[2] );
        if( length == 3 ) return result;
        list strlist(_ob[3]);
        while( len( strlist ) )
          result->type.push_back( extract<std::string>(strlist.pop(0)) );
        if( length == 4 ) return result;
        result->site = extract< types::t_real >( _ob[4] );
        return result;
      }
      catch( std::exception &_e )
      {
        delete result;
        std::ostringstream sstr;
        sstr << "Object cannot be converted to an atom: \n"
             << _e.what() << "\n";
        throw std::runtime_error( sstr.str() );
      }
      catch( ... )
      {
        delete result;
        throw std::runtime_error( "Could not convert object to Atom." );
      }
      return NULL;
    }
    types::t_real toReal(std::string _str )
    { 
      if( not Crystal::Structure::lattice )
        throw std::runtime_error( "Could not convert atom type.\n" ); 
      if( _str.compare( Crystal::Structure::lattice->sites[0].type[0] ) == 0 )
        return types::t_real(-1);
      else if( _str.compare( Crystal::Structure::lattice->sites[0].type[1] ) == 0 )
        return types::t_real(1);
      else
        throw std::runtime_error( "Requested Atomic type is not within lattice" );
    }
    std::string toType( types::t_real _r )
    { 
      if( not Crystal::Structure::lattice )
        throw std::runtime_error( "Could not convert atom type.\n" ); 
      if( Crystal::Structure::lattice->sites.size() != 1)
        throw std::runtime_error( "Lattice cannot have more than one atomic site.\n" ); 
      if( Crystal::Structure::lattice->sites[0].type.size() != 2)
        throw std::runtime_error( "Lattice must have two atomic types.\n" ); 
      if(     Fuzzy::neq( _r, types::t_real(1) ) 
          and Fuzzy::neq( _r, types::t_real(-1) ) ) return "";
      return Fuzzy::neq( _r, types::t_real(1) ) ? 
                Crystal::Structure::lattice->sites[0].type[0]:
                Crystal::Structure::lattice->sites[0].type[1] ;
    }

  }
} // namespace LaDa
