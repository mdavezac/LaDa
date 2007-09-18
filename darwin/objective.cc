//
//  Version: $Id$
//
#include  "objective.h"

namespace GA
{
  bool Fitness :: Load( const TiXmlElement & _node )
  {
    if ( not _node.Attribute( "fitness" ) )
    {
      is_valid = false;
      return false; 
    }
    double d; _node.Attribute( "fitness", &d );
    quantity = (t_Quantity) d;
    is_valid = true;
    return true;
  }
  bool Fitness :: Save( TiXmlElement & _node ) const
  {
    if ( not is_valid )
    {
      std::cerr << "Trying to Save invalid fitness" << std::endl;
      return false;
    }
    double d = (double) quantity;
    _node.SetDoubleAttribute("fitness", d);
    return true;
  }
}


