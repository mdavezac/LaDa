//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/error.hpp>

#include <opt/types.h>
#include <opt/debug.h>

#include "ce.hpp"
#include "ce.impl.hpp"

namespace LaDa
{
  namespace Python
  {
    template<class T_HARMONIC>
      void expose_ce_functional( const std::string &_name, const std::string &_docstring );

    template<class T_HARMONIC>
      void load_builder( CE::Builder<T_HARMONIC>& _functional, const std::string &_filename )
      {
        TiXmlDocument doc( _filename ); 
        TiXmlHandle docHandle( &doc ); 
      
        if( not doc.LoadFile() ) 
        {
          PyErr_SetString
          (
            PyExc_IOError, 
            ("Could not open/parse " + _filename + ": " + doc.ErrorDesc() + "\n" ).c_str() 
          );
          bp::throw_error_already_set();
          return;
        }
        const TiXmlElement* parent = docHandle.FirstChild("Job").Element();
        if( not parent )
        {
          PyErr_SetString
          (
            PyExc_IOError, 
            ("Could not find <Job> </Job> tags in " + _filename + "\n").c_str() 
          );
          bp::throw_error_already_set();
          return;
        }
      
        try
        {
          if( not _functional.Load( *parent ) )
          {
            PyErr_SetString
            ( 
              PyExc_IOError, 
              ("Could not load ce functional from " + _filename + "\n").c_str() 
            );
            bp::throw_error_already_set();
            return;
          }
          _functional.add_equivalent_clusters();
        }
        catch( std::exception &_e )
        {
          PyErr_SetString
          (
            PyExc_IOError, 
            ("Could not create ce functional from input-file " + _filename + "\n").c_str() 
          );
          bp::throw_error_already_set();
        }
      }

    template< class T_HARMONIC >
      std::pair
      <
        typename CE::Builder<T_HARMONIC>::t_Chemical*,
        typename CE::Builder<T_HARMONIC>::t_CS*
      > create( const CE::Builder<T_HARMONIC> &_functional, const Crystal::Structure &_str )
      {
        if( ptr_str.atoms.size() ) 
        {
          PyErr_SetString( PyExc_RuntimeError, "Structure is empty.\n" );
          bp::throw_error_already_set();
        }
        try
        {
          if ( _str.kvecs.size() == 0 )
          {
            std::auto_ptr<Crystal::Structure> str( new Crystal::Structure( _str ) );
            str.find_k_vectors();
            call( _functional, *str );
          }
        }
        catch(...) { bp::throw_error_already_set(); return }

        try { return _functional.generate_functional( _str ); }
        catch( std::exception &_e )
        {
          PyErr_SetString( PyExc_IOError, "Could not evaluate CE functional.\n" ); 
          bp::throw_error_already_set();
        }
      }

    template< class T_HARMONIC >
      types::t_real create( const CE::Builder<T_HARMONIC> &_functional,
                            const Crystal::Structure &_str )
      {
        typedef std::pair 
                <
                  typename CE::Builder<T_HARMONIC>::t_Chemical*,
                  typename CE::Builder<T_HARMONIC>::t_CS*
                > t_Pair;
      }

    namespace XML
    {
      template<> std::string nodename<t_CubicCS>()    { return "CS"; }
      template<> std::string nodename<t_TetraCS>()    { return "CS"; }
      
      template<> void do_specialcode< t_CubicBuilder >( t_CubicBuilder &_type )
        { _type.add_equivalent_clusters(); }
      template<> void do_specialcode< t_TetraBuilder >( t_TetraBuilder &_type )
        { _type.add_equivalent_clusters(); }
      
      template<> bool doloadcode<t_CubicCS>( t_CubicCS &_type, TiXmlElement *_parent )
        { return _type.Load_Harmonics( *_parent ); }
      template<> bool doloadcode<t_TetraCS>( t_TetraCS &_type, TiXmlElement *_parent )
        { return _type.Load_Harmonics( *_parent ); }
        
      template<> TiXmlElement *findnode<t_CubicCS>( TiXmlHandle &_doc )
      {
        TiXmlElement *result = _doc.FirstChild("Job")
                                   .FirstChild(nodename<t_CubicCS>()).Element();
        if( result ) return result;
        
        result = _doc.FirstChild("Job")
                     .FirstChild("Functional").Element();
        for(; result; result = result->NextSiblingElement("Functional") )
        {
          if( not result->Attribute("type") ) continue;
        
          std::string strg = result->Attribute("type");
          if( strg.compare("CE") == 0 ) break;
        }   
        if( not result ) return NULL;
        result = result->FirstChildElement(nodename<t_CubicCS>());
        return result;
      }
      template<> TiXmlElement *findnode<t_TetraCS>( TiXmlHandle &_doc )
        { return findnode<t_CubicCS>( _doc ); }
      template<> TiXmlElement *findnode<t_CubicBuilder>( TiXmlHandle &_doc )
        { return _doc.FirstChild("Job").Element(); }
      template<> TiXmlElement *findnode<t_TetraBuilder>( TiXmlHandle &_doc )
        { return findnode<t_CubicBuilder>( _doc ); }
    } // end of XML namespace

    void expose_ce()
    {
      using namespace boost::python;

      typedef CE::Builder< CE::ConstituentStrain::Harmonic::Cubic > :: t_Chemical t_Chemical;
      class_< t_Chemical >( "Chemical" )
        .def( init< t_Chemical >() )
        .def( "evaluate", &t_Chemical::evaluate );

      details::ExposeHarmonicRelated< CE::ConstituentStrain::Harmonic::Cubic >();
      details::ExposeHarmonicRelated< CE::ConstituentStrain::Harmonic::Tetragonal >();
    }
  } // end of PythonLaDa namespace
} // namespace LaDa
