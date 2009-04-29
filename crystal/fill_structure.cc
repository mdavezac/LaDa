//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <functional>

#include <opt/ndim_iterator.h>
#include <opt/smith_normal_form.h>
#include <opt/debug.h>
#include <atat/misc.h>

#include "lattice.h"
#include "structure.h"
#include "fill_structure.h"



namespace LaDa
{

  namespace Crystal 
  {
    bool fill_structure( Crystal::Structure &_structure )
    {
      __ASSERT( _structure.lattice == NULL, "Lattice not set.\n" )
      const atat::rMatrix3d inv_lat( !_structure.lattice->cell );
      const atat::rMatrix3d inv_lat_cell( inv_lat * _structure.cell );
      atat::iMatrix3d int_cell;
      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
        {
          const types::t_real d( std::floor( inv_lat_cell(i,j) + 0.5  ) );
          int_cell(i,j) = types::t_int( d );
          __DOASSERT
          ( 
            std::abs( types::t_real( int_cell(i,j) ) - inv_lat_cell(i,j) ) > 0.01,
               "Input structure is not supercell of the lattice: " 
            << int_cell(i,j) << " != " << inv_lat_cell(i,j) << "\n"
          )
        }
      atat::iMatrix3d left, right, smith;
      opt::smith_normal_form( smith, left, int_cell, right );

      atat::rMatrix3d rleft;
      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
          rleft(i,j) = types::t_real( left(i,j) );
      const atat::rMatrix3d factor
      ( 
        ( !_structure.cell ) * _structure.lattice->cell * (!rleft) 
      );
      opt::Ndim_Iterator< types::t_int, std::less<types::t_int> > iterator;
      iterator.add( 0, smith(0,0) );
      iterator.add( 0, smith(1,1) );
      iterator.add( 0, smith(2,2) );
      
      do
      {
        // in "quotient" space
        const atat::rVector3d vec1
        ( 
          ( types::t_real ) iterator.access(0),
          ( types::t_real ) iterator.access(1),
          ( types::t_real ) iterator.access(2)
        );
        // in supercell fractional
        const atat::rVector3d vec2( factor * vec1 );
        // in supercell fractional and in supercell parallelogram
        const atat::rVector3d vec3
        (
          vec2(0) - std::floor( vec2(0) ),
          vec2(1) - std::floor( vec2(1) ),
          vec2(2) - std::floor( vec2(2) )
        );
        // in cartesian
        const atat::rVector3d vec( _structure.cell * vec3 );

        // adds all lattice sites.
        typedef Crystal::Lattice::t_Site t_Site;
        foreach( const t_Site &site, _structure.lattice->sites ) 
        {
          Structure::t_Atom atom;
          atom.pos = vec + site.pos;
          atom.site = site.site;
          atom.freeze = site.freeze;
          atom.type = -1e0;
          _structure.atoms.push_back(atom);
        }
      } while( ++iterator );

      return true;
    }


//   bool fill_structure2( Crystal::Structure &_structure )
//   {
//     if( not _structure.lattice ) return false;
//     
//     Structure :: t_Atoms atoms;
//     atat::rVector3d vec;
//     atat::rMatrix3d &cell = _structure.cell; 
//     Crystal::Lattice &lattice = *_structure.lattice; 
//
//     // Construct the transition matrix from the lattice basis to this cell-shape basis
//     atat::rMatrix3d M = (!cell) * lattice.cell;
//     // The next few operations should tell us the maximum range of the structure
//     // cell int terms of the lattice cell.
//     atat::rVector3d r = (!lattice.cell) * ( cell * atat::rVector3d(1,1,1) );
//     types::t_int range =   (types::t_int) std::ceil( atat::det( cell )
//                          / atat::det(lattice.cell) );
//     
//     // now that we have the range, we can fully explore the region
//     // sets up the n-dimensional iterators
//     opt::Ndim_Iterator< types::t_int, std::less_equal<types::t_int> > global_iterator;
//     global_iterator.add( -range, range);
//     global_iterator.add( -range, range);
//     global_iterator.add( -range, range);
//
//     do
//     {
//       // creates vector in lattice basis
//       vec[0] =  (types::t_real) global_iterator.access(0);
//       vec[1] =  (types::t_real) global_iterator.access(1);
//       vec[2] =  (types::t_real) global_iterator.access(2);
//       // Transforms vec to fractional coordinates in current cell-shape
//       vec = M * vec;
//       // if any of the coordinates is >= 1, then this is a periodic image
//       if (    Fuzzy::geq( vec(0), 1.0 ) 
//            or Fuzzy::geq( vec(1), 1.0 ) 
//            or Fuzzy::geq( vec(2), 1.0 ) ) continue;
//       // if any of the coordinates is < 0, then this is a periodic image
//       if (    Fuzzy::le( vec(0), 0.0 ) 
//            or Fuzzy::le( vec(1), 0.0 ) 
//            or Fuzzy::le( vec(2), 0.0 ) ) continue;
//
//
//       // Goes back to lattice basis
//       vec[0] =  (types::t_real) global_iterator.access(0);
//       vec[1] =  (types::t_real) global_iterator.access(1);
//       vec[2] =  (types::t_real) global_iterator.access(2);
//       // And then to cartesian
//       vec = lattice.cell * vec;
//
//       atoms.push_back( Crystal::Structure::t_Atom(vec,0) );
//       
//     } while( ++global_iterator );
//
//     // Finally, we copy the kvector positions as atoms, and the related sites
//     Structure::t_Atoms::const_iterator i_vec = atoms.begin();
//     Structure::t_Atoms::const_iterator i_vec_end = atoms.end();
//     _structure.atoms.clear();
//     _structure.atoms.reserve( _structure.lattice->sites.size() * atoms.size() );
//     bool only_one_site = _structure.lattice->sites.size() == 1;
//     Lattice :: t_Sites :: const_iterator i_site;
//     Lattice :: t_Sites :: const_iterator i_site_end = _structure.lattice->sites.end();
//     for(; i_vec != i_vec_end; ++i_vec )
//       for( i_site = _structure.lattice->sites.begin(); i_site != i_site_end; ++i_site )
//       {
//         Structure::t_Atom atom;
//         atom.pos = i_vec->pos + i_site->pos;
//         atom.site = i_site->site;
//         atom.freeze = i_site->freeze;
//         atom.type = -1e0;
//         _structure.atoms.push_back(atom);
//       }
//     return true;
//   }
   

  } // namespace Crystal

} // namespace LaDa
