//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <print/stdout.h>
#include <print/manip.h>
#include <print/xmg.h>

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/string.hpp>

#include "bandgap_stubs.h"

namespace LaDa
{
  namespace GA
  {
    namespace Keepers 
    {
      bool BandGap :: Load ( const TiXmlElement &_node )
      {
        double d;

        if ( not _node.Attribute("cbm", &d ) ) goto errorout;
        cbm = types::t_real(d);
        if ( not _node.Attribute("vbm", &d ) ) goto errorout;
        vbm = types::t_real(d);

        return true;

        errorout:
          std::cerr << "Could not Load GA::Keepers::BandGap" << std::endl;
          return false;
      }
      bool BandGap :: Save( TiXmlElement &_node ) const
      {
        _node.SetDoubleAttribute("vbm", vbm );
        _node.SetDoubleAttribute("cbm", cbm );
        if( not Vff::Keeper::Save( _node ) ) return false;

        return true;
      }
    } // namespace Keepers
  } // namespace GA
}


