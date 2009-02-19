//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <opt/types.h>
#include <opt/debug.h>

#include "ce.hpp"
#include "ce.impl.hpp"

namespace LaDa
{
  namespace Python
  {
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

      class_< std::vector<types::t_real> >( "details_variables" )
        .def( vector_indexing_suite< std::vector<types::t_real> >() );
        
      details::ExposeHarmonicRelated< CE::ConstituentStrain::Harmonic::Cubic >();
      details::ExposeHarmonicRelated< CE::ConstituentStrain::Harmonic::Tetragonal >();
    }
  } // end of PythonLaDa namespace
} // namespace LaDa
