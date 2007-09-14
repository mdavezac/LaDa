//
//  Version: $Id: pescan.cc 275 2007-09-07 17:42:43Z davezac $
//
#include <functional>
#include <algorithm>
#include <ext/algorithm>
#include <fstream>
#ifndef __PGI
  #include<ext/functional>
  using __gnu_cxx::compose1;
#else
  #include<functional>
  using std::compose1;
#endif

#include <print/stdout.h>
#include <print/manip.h>
#include <lamarck/structure.h>
#include <lamarck/atom.h>
#include <opt/va_minimizer.h>

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

  void Evaluator::init( t_Individual &_indiv )
  {
    t_Base :: init( _indiv );
    // sets structure to this object 
    structure << *current_object;
  }
  void Evaluator::evaluate()
  {
    // relax structure
    vff();
    // Load relaxed structure into pescan
    pescan << vff; 
    // get band gap
    pescan( *current_object );

    // set quantity
    current_individual->quantities() = current_object->cbm - current_object->vbm;
  }

} // namespace pescan



