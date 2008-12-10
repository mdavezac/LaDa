//
//  Version: $Id$
//

namespace LaDa
{
  namespace GA
  {
    namespace PureLayers
    {

#     if defined( EVALHEAD ) || defined(INEVAL)
#       error "Macros with same names."
#     endif
#     define EVALHEAD  Evaluator<T_INDIVIDUAL, T_TRANSLATE, T_ASSIGN> 
#     define INEVAL( var ) \
         template< class T_INDIVIDUAL, \
                   template<class> class T_TRANSLATE, \
                   template<class,class> class T_ASSIGN >  var EVALHEAD

      INEVAL( bool ) :: initialize( t_Individual &_indiv )
      {
        foreach( Crystal::Structure::t_Atom &atom, structure.atoms )
          if ( i_atom->freeze & Crystal::Structure::t_Atom::FREEZE_T ) 
            atom.type = eo::rng.flip() ? 1e0: -1e0;
        translate( structure, _indiv.Object() );
        concentration( _indiv.Object() );
        _indiv.invalidate();
        return true;
      }
      INEVAL(bool) :: Load( const TiXmlElement &_node )
      {
        if( not t_Base :: Load( _node ) ) return false;
        return concentration.Load( _node );
      }


#     undef EVALHEAD
#     undef INEVAL

    } // namespace Layered


 }
} // namespace LaDa
 

#include "evaluator.impl.h"

#endif // _LAYERED_H_
