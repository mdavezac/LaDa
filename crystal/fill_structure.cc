//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/ndim_iterator.h>
#include <opt/debug.h>

#include <atat/misc.h>

#include "lattice.h"
#include "structure.h"
#include "fill_structure.h"



namespace LaDa
{

  namespace Crystal 
  {

    bool fill_structure( Crystal::Structure &_str )
    {
      if( not _str.lattice ) return false;
      
      Structure :: t_Atoms atoms;
      atat::rVector3d vec;
      atat::rMatrix3d &cell = _str.cell; 
      Crystal::Lattice &lattice = *_str.lattice; 

      // Construct the transition matrix from the lattice basis to this cell-shape basis
      atat::rMatrix3d M = (!cell) * lattice.cell;
      // The next few operations should tell us the maximum range of the structure
      // cell int terms of the lattice cell.
      atat::rVector3d r = (!lattice.cell) * ( cell * atat::rVector3d(1,1,1) );
      types::t_int range =   (types::t_int) std::ceil( atat::det( cell )
                           / atat::det(lattice.cell) );
      
      // now that we have the range, we can fully explore the region
      // sets up the n-dimensional iterators
      opt::Ndim_Iterator< types::t_int, std::less_equal<types::t_int> > global_iterator;
      global_iterator.add( -range, range);
      global_iterator.add( -range, range);
      global_iterator.add( -range, range);

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

        atoms.push_back( Crystal::Structure::t_Atom(vec,0) );
        
      } while( ++global_iterator );

      // Finally, we copy the kvector positions as atoms, and the related sites
      Structure::t_Atoms::const_iterator i_vec = atoms.begin();
      Structure::t_Atoms::const_iterator i_vec_end = atoms.end();
      _str.atoms.clear();
      _str.atoms.reserve( _str.lattice->sites.size() * atoms.size() );
      bool only_one_site = _str.lattice->sites.size() == 1;
      Lattice :: t_Sites :: const_iterator i_site;
      Lattice :: t_Sites :: const_iterator i_site_end = _str.lattice->sites.end();
      for(; i_vec != i_vec_end; ++i_vec )
        for( i_site = _str.lattice->sites.begin(); i_site != i_site_end; ++i_site )
        {
          Structure::t_Atom atom;
          atom.pos = i_vec->pos + i_site->pos;
          atom.site = i_site->site;
          atom.freeze = i_site->freeze;
          atom.type = -1e0;
          _str.atoms.push_back(atom);
        }
      return true;
    }
   

  } // namespace Crystal

} // namespace LaDa
