//
//  Version: $Id$
//
#include "bandgap.h"

namespace BandGap
{
  bool Evaluator :: Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type )
  {
    t_Object &object = _indiv.Object();
    if ( not object.Load( _node ) ) return false;
    _indiv.quantities() = object.cbm - object.vbm; 

    return t_Base::Load( _indiv, _node, _type );
  }
  bool Evaluator :: Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const
  { 
    return _indiv.Object().Save(_node) and t_Base::Save( _indiv, _node, _type );
  }
  bool Evaluator :: Load( const TiXmlElement &_node )
  {
    __DOASSERT( not t_Base::Load( _node ),
                " Could not load TwoSites::Evaluator<Object> input.\n" )
    __DOASSERT( not bandgap.Load( _node ),
                " Could not load bandgap interface from input.\n" )

    return true;
  }

  void Evaluator::evaluate()
  {
    Crystal::Structure copy_structure = structure;
    concentration.get( *current_object );
    current_object->x = concentration.x;
    current_object->y = concentration.y;
    // get band gap
#ifdef _NOLAUNCH
    typedef t_Individual :: t_IndivTraits :: t_FourierRtoK t_Fourier;
    t_Fourier( structure.atoms.begin(), structure.atoms.end(),
               structure.k_vecs.begin(), structure.k_vecs.end() );
#endif
    bandgap( *current_object );

    // set quantity
    current_individual->quantities() = current_object->cbm - current_object->vbm;
    structure = copy_structure;
  }

} // namespace BandGap



