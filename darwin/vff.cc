//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "vff.h"

namespace Vff
{
  bool Keeper :: Load( const TiXmlElement &_node )
  {
    const TiXmlElement *vffxml = _node.FirstChildElement( "VffResult" );
    atat::rVector3d vec;
    double d;

    if ( not vffxml ) goto errorout;
    
    if ( not vffxml->Attribute("energy", &d ) ) goto errorout;
    energy = types::t_real(d);
    if ( not vffxml->Attribute("xx",     &d ) ) goto errorout;
    vec(0) = types::t_real(d);
    if ( not vffxml->Attribute("xy",     &d ) ) goto errorout;
    vec(1) = types::t_real(d);
    if ( not vffxml->Attribute("xz",     &d ) ) goto errorout;
    vec(2) = types::t_real(d);
    stress.set_column(0,vec);
    if ( not vffxml->Attribute("yx",     &d ) ) goto errorout;
    vec(0) = types::t_real(d);
    if ( not vffxml->Attribute("yy",     &d ) ) goto errorout;
    vec(1) = types::t_real(d);
    if ( not vffxml->Attribute("yz",     &d ) ) goto errorout;
    vec(2) = types::t_real(d);
    stress.set_column(1,vec);
    if ( not vffxml->Attribute("zx",     &d ) ) goto errorout;
    vec(0) = types::t_real(d);
    if ( not vffxml->Attribute("zy",     &d ) ) goto errorout;
    vec(1) = types::t_real(d);
    if ( not vffxml->Attribute("zz",     &d ) ) goto errorout;
    vec(2) = types::t_real(d);
    stress.set_column(2,vec);

    return true;
errorout:
    std::cerr << "Could not Load Vff::Keeper" << std::endl;
    return false;
  }
  bool Keeper :: Save( TiXmlElement &_node ) const
  {
    TiXmlElement *vffxml = new TiXmlElement( "VffResult" );
    if ( not vffxml )
    {
      std::cerr << "Could not Save Vff::Keeper";
      return false;
    }
    vffxml->SetDoubleAttribute("energy", energy );
    vffxml->SetDoubleAttribute("xx", stress.get_column(0)(0));
    vffxml->SetDoubleAttribute("xy", stress.get_column(0)(1));
    vffxml->SetDoubleAttribute("xz", stress.get_column(0)(2));
    vffxml->SetDoubleAttribute("yx", stress.get_column(1)(0));
    vffxml->SetDoubleAttribute("yy", stress.get_column(1)(1));
    vffxml->SetDoubleAttribute("yz", stress.get_column(1)(2));
    vffxml->SetDoubleAttribute("zx", stress.get_column(2)(0));
    vffxml->SetDoubleAttribute("zy", stress.get_column(2)(1));
    vffxml->SetDoubleAttribute("zz", stress.get_column(2)(2));

    _node.LinkEndChild( vffxml );

    return true;
  }

} // namespace pescan



