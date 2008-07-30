//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gsl_simplex.h"

namespace Minimizer {



  void Simplex :: set_parameters( types::t_unsigned _itermax,
                                  types::t_real _tol, 
                                  types::t_real _stepsize ) 
  {
    tolerance     = _tol;
    stepsize      = _stepsize;
    itermax       = _itermax;
  }
            

  const TiXmlElement* Simplex :: find_node( const TiXmlElement &_node )
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
      if ( not parent->Attribute( "type" )  ) break;
      str = parent->Attribute("type");
      if ( str.compare("simplex" ) == 0 ) break;
      parent = parent->NextSiblingElement("Minimizer");
    }
    if ( parent ) return parent;
    
    return NULL;
  }
  
  bool Simplex :: Load_( const TiXmlElement &_node )
  {
    double d; int n;
    _node.Attribute( "itermax", &n );
    itermax = (n > 0) ? abs(n) : 10000;
    _node.Attribute( "tolerance", &d );
    tolerance = d ? types::t_real(d) : types::tolerance;
    _node.Attribute( "stepsize", &d );
    stepsize = d ? types::t_real(d) : 0.1;
    if( _node.Attribute("verbose") ) verbose = true;
  
    return true;
  }
  
  bool Simplex :: Load( const TiXmlElement &_node )
  {
    const TiXmlElement* parent = find_node( _node );
    if( parent ) return Load_(*parent);
    std::cerr << "Could not find an <Minimizer type=\"simplex\"> tag in input file" 
              << std::endl;
    return false;
  }
}

