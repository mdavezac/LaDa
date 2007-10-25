//
//  Version: $Id$
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

#include "molecularity.h"
#include "two_sites.h"

  // "Molecularity" means that we try and simulate growth conditions...
  // The input should be a 1x1xn supercell, with n either 100, 110 or 111 
  // each unit-cell is a "molecule", say InAs or GaSb, but not a mix of the two.


namespace Molecularity
{
  bool Evaluator :: Load( t_Individual &_indiv,
                          const TiXmlElement &_node, 
                          bool _type )
  {
    t_Object &object = _indiv.Object();
    if ( not object.Load( _node ) ) return false;

    // set quantity
    object_to_quantities( _indiv );

    return t_Base::Load( _indiv, _node, _type );
  }
  
  bool Evaluator :: Load( const TiXmlElement &_node )
  {
    if ( not t_Base::Load( _node ) ) return false;

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

} // namespace Molecularity



