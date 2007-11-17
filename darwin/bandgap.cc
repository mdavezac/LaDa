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
    if ( not t_Base::Load( _node ) )
    {
      std::cerr << " Could not load TwoSites::Evaluator<Object> input!! " << std::endl; 
      return false;
    }
    if ( not vff.Load( _node ) )
    {
      std::cerr << " Could not load vff input!! " << std::endl; 
      return false;
    }
    if ( not pescan.Load( _node ) )
    {
      std::cerr << " Could not load pescan interface from input!! " << std::endl; 
      return false;
    }

    return true;
  }

  void Evaluator::evaluate()
  {
    concentration.get( *current_object );
    current_object->x = concentration.x;
    current_object->y = concentration.y;
    // relax structure
    vff();
    // Load relaxed structure into pescan
    pescan << vff; 
    // get band gap
#ifdef _NOLAUNCH
    typedef t_Individual :: t_IndivTraits :: t_FourierRtoK t_Fourier;
    t_Fourier( structure.atoms.begin(), structure.atoms.end(),
               structure.k_vecs.begin(), structure.k_vecs.end() );
#endif
    pescan( *current_object );

    // set quantity
    current_individual->quantities() = current_object->cbm - current_object->vbm;
  }

} // namespace pescan



