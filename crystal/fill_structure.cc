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
      const atat::rMatrix3d inv_lat_cell( ( !_structure.lattice->cell ) * _structure.cell );
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
      Structure::t_Atom atom; atom.type = -1e0;
      for( size_t i(0); i < smith(0,0); ++i )
        for( size_t j(0); j < smith(1,1); ++j )
          for( size_t k(0); k < smith(2,2); ++k )
          {
            // in supercell fractional
            const atat::rVector3d vec1( factor * atat::rVector3d(i,j,k) );
            // in supercell fractional and in supercell parallelogram
            const atat::rVector3d vec2
            (
              vec1(0) - std::floor( vec1(0) ),
              vec1(1) - std::floor( vec1(1) ),
              vec1(2) - std::floor( vec1(2) )
            );
            // in cartesian
            const atat::rVector3d vec( _structure.cell * vec1 );
          
            // adds all lattice sites.
            typedef Crystal::Lattice::t_Site t_Site;
            foreach( const t_Site &site, _structure.lattice->sites ) 
            {
              atom.pos = vec + site.pos;
              atom.site = site.site;
              atom.freeze = site.freeze;
              _structure.atoms.push_back(atom);
            }
          }

      return true;
    }


   

  } // namespace Crystal

} // namespace LaDa
