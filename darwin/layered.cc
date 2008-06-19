//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>
#include <opt/ndim_iterator.h>

#include "layered.h"

void FillStructure( Crystal::Structure &_str )
{
  if( not _str.lattice ) return;
  
  atat::rVector3d vec;
  atat::rMatrix3d &cell = _str.cell; 
  Crystal::Lattice &lattice = *_str.lattice; 

  // Construct the transition matrix from the lattice basis to this cell-shape basis
  atat::rMatrix3d M = (!cell) * lattice.cell;
  // The next few operations should tell us the maximum range of the structure
  // cell int terms of the lattice cell.
  atat::rVector3d r = (!lattice.cell) * ( cell * atat::rVector3d(1,1,1) );
  atat::iVector3d range( (types::t_int) std::ceil( std::abs(r(0) ) ), 
                         (types::t_int) std::ceil( std::abs(r(1) ) ), 
                         (types::t_int) std::ceil( std::abs(r(2) ) ) );
  
  // now that we have the range, we can fully explore the region
  // sets up the n-dimensional iterators
  opt::Ndim_Iterator< types::t_int, std::less_equal<types::t_int> > global_iterator;
  global_iterator.add( -range[0], range[0]);
  global_iterator.add( -range[1], range[1]);
  global_iterator.add( -range[2], range[2]);
//   std::cout << (!lattice.cell) * cell << std::endl
//             << "range: " << r << "  " << range << std::endl;

  
  _str.atoms.clear();
  do
  {
    // creates vector in lattice basis
    vec[0] =  (types::t_real) global_iterator.access(0);
    vec[1] =  (types::t_real) global_iterator.access(1);
    vec[2] =  (types::t_real) global_iterator.access(2);
    // Transforms vec to fractional coordinates in current cell-shape
    vec = M * vec;
    // if any of the coordinates is >= 1, then this is a periodic image
    if (    Fuzzy::geq( vec(0), 1.0 ) 
         or Fuzzy::geq( vec(1), 1.0 ) 
         or Fuzzy::geq( vec(2), 1.0 ) ) continue;
    // if any of the coordinates is < 0, then this is a periodic image
    if (    Fuzzy::le( vec(0), 0.0 ) 
         or Fuzzy::le( vec(1), 0.0 ) 
         or Fuzzy::le( vec(2), 0.0 ) ) continue;


    // Goes back to lattice basis
    vec[0] =  (types::t_real) global_iterator.access(0);
    vec[1] =  (types::t_real) global_iterator.access(1);
    vec[2] =  (types::t_real) global_iterator.access(2);
    // And then to cartesian
    vec = lattice.cell * vec;

    _str.atoms.push_back( Crystal::Structure::t_Atom(vec,0) );
    
  } while( ++global_iterator );

}
