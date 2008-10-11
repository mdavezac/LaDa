//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <functional>
#include <algorithm>

#include <crystal/structure.h>

#include "groundstate.h"
#include "functors.h"
#include "concentration.h"

namespace GroundState
{
  bool Evaluator :: Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type )
  {
    double d;
    if ( not _node.Attribute("CE", &d ) ) goto errorout;
    _indiv.quantities() = (types::t_real) d;

    return t_Base::Load( _indiv, _node, _type );
errorout:
    std::cerr << "Could not Load CE::Object" << std::endl;
    return false;
  }
  bool Evaluator :: Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const
  { 
    double d = (double) _indiv.const_quantities();
    _node.SetDoubleAttribute("CE", d );
    return t_Base::Save( _indiv, _node, GA::LOADSAVE_SHORT );
  }

  bool Evaluator :: Load( const TiXmlElement &_node )
  {
    if( not t_Base::Load( _node ) ) return false;
    
    if ( not ce.Load( _node ) ) return false;
    
    return true;
  }

} // namespace GroundState


