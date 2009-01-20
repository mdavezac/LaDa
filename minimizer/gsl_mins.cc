//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/lexical_cast.hpp>

#include "gsl_mins.h"

namespace LaDa
{
  namespace Minimizer {

    LADA_REGISTER_MINIMIZER_VARIANT_SOURCE( Gsl, "GSL" )


    void Gsl :: set_parameters( t_gsl_minimizer_type _type, 
                                types::t_unsigned _itermax,
                                types::t_real _tol, 
                                types::t_real _linetol, 
                                types::t_real _linestep ) 
    {
      tolerance     = _tol;
      linetolerance = _linetol;
      linestep      = _linestep;
      itermax       = _itermax;
      type          = _type;
    }
              

    const TiXmlElement* Gsl :: find_node( const TiXmlElement &_node )
    {
      const TiXmlElement *parent;
      std::string str;
    
      // This whole section tries to find a <Functional type="vff"> tag
      // in _node or its child
      str = _node.Value();
      if ( str.compare("Minimizer" ) != 0 )
        parent = _node.FirstChildElement("Minimizer");
      else parent = &_node;
    
      
      while (parent)
      {
        str = "";
        type = GSL_NONE;
        if ( parent->Attribute( "type" )  )
          str = parent->Attribute("type");
        if ( str.compare("gsl_bfgs2" ) == 0 )     type = GSL_BFGS2;
        else if ( str.compare("gsl_fr" ) == 0 )   type = GSL_FR;
        else if ( str.compare("gsl_pr" ) == 0 )   type = GSL_PR;
        else if ( str.compare("gsl_bfgs" ) == 0 ) type = GSL_BFGS;
        else if ( str.compare("gsl_sd" ) == 0 )   type = GSL_SD;
        if ( type != GSL_NONE ) break;
        parent = parent->NextSiblingElement("Minimizer");
      }
      if ( parent ) return parent;
      
      return NULL;
    }
    
    bool Gsl :: Load_( const TiXmlElement &_node )
    {
      if( _node.Attribute( "itermax" ) )
        itermax = boost::lexical_cast<types::t_unsigned>( _node.Attribute("itermax") );
      if( _node.Attribute( "tolerance" ) )
        tolerance = boost::lexical_cast<types::t_real>( _node.Attribute("tolerance") );
      if( _node.Attribute( "linetolerance" ) )
        linetolerance = boost::lexical_cast<types::t_real>( _node.Attribute("linetolerance") );
      if( _node.Attribute( "linestep" ) )
        linestep = boost::lexical_cast<types::t_real>( _node.Attribute("linestep") );
      if( _node.Attribute("verbose") ) 
      {
        const std::string value( _node.Attribute("verbose") );
        if( value == "true" or value == "TRUE" or value == "T" or value == "t" ) 
          verbose = true;
        else if( value == "false" or value == "FALSE" or value == "F" or value == "f" ) 
          verbose = false;
        else verbose = boost::lexical_cast<bool>( _node.Attribute("verbose") );
      }
      return true;
    }
    
    bool Gsl :: Load( const TiXmlElement &_node )
    {
      const TiXmlElement* parent = find_node( _node );
      if( parent ) return Load_(*parent);
      return false;
    }
  }
} // namespace LaDa
