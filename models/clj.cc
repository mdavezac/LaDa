//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "clj.h"

namespace LaDa
{
  namespace Models
  {
    bool Clj :: Load( const TiXmlElement& _node )
    {
      const TiXmlElement* const parent = opt::find_node( _node, "Functional", "Clj" );
      if( not parent ) return false;
      if( not Ewald::Load( *parent ) ) return false;
      if( not LennardJones :: Load( *parent ) ) return false;
      return true;
    }

    Clj :: t_Return Clj :: energy(const t_Arg& _in, t_Arg& _out) const
    {
      _out.cell.zero();
      foreach( t_Arg :: t_Atom &atom, _out.atoms )
        atom.pos = atat::rVector3d(0,0,0);

      _out.energy = Ewald::energy( _in, _out ) + LennardJones::energy( _in, _out );
      return _out.energy;
    }
  } // namespace CLJ
} // namespace LaDa
