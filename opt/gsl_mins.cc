//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gsl_mins.h"

namespace Minimizer {



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
    double d; int n;
    _node.Attribute( "itermax", &n );
    itermax = (n > 0) ? abs(n) : 10000;
    _node.Attribute( "tolerance", &d );
    tolerance = d ? types::t_real(d) : types::tolerance;
    _node.Attribute( "linetolerance", &d );
    linetolerance = d ? types::t_real(d) : 0.01;
    _node.Attribute( "linestep", &d );
    linestep = d ? types::t_real(d) : 0.1;
    if( _node.Attribute("verbose") ) verbose = true;
  
    return true;
  }
  
  bool Gsl :: Load( const TiXmlElement &_node )
  {
    const TiXmlElement* parent = find_node( _node );
    if( parent ) return Load_(*parent);
    std::cerr << "Could not find an <Minimizer type=\"some gsl\"> tag in input file" 
              << std::endl;
    return false;
  }
}

