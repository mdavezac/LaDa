//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/debug.h>
#include <opt/types.h>

#include <Eigen/LU>

#include "lattice.h"
#include "structure.h"
#include "epi_structure.h"
#include "fill_structure.h"



namespace LaDa
{

  namespace Crystal 
  {
    bool create_epitaxial_structure( Structure& _structure,
                                     math::rVector3d &_direction,
                                     math::iVector3d &_extent )
    {

      // Then constructs unit cell
      math::rMatrix3d &cell( _structure.cell ); 
      math::rMatrix3d diagonal( _extent.cast<types::t_real>().asDiagonal() );
      Lattice *lattice( _structure.lattice );
      cell = lattice->cell;
      cell.col(0) = lattice->cell * _direction; 
      cell = cell * diagonal;

      // Checks that cell is not singular
      if ( math::is_zero( std::abs(cell.determinant()) ) )
        cell.col(1) = lattice->cell.col( 0 );
      if ( math::is_zero( std::abs(cell.determinant()) ) )
        cell.col(2) = lattice->cell.col( 1 );
      if ( math::is_zero( std::abs(cell.determinant()) ) )
      {
        std::cerr << "Could not construct unit-cell\n" << cell << std::endl;
        return false;
      }
      if( math::is_zero( cell.determinant() ) )
      {
        math::rVector3d swap = cell.col(1);
        cell.col(1) = cell.col(2);
        cell.col(2) = swap;
      }

      // Makes sure the triad is direct
      if ( cell.determinant() < 0 )
      {
        math::rVector3d d = cell.col(2);
        cell.col(2) = cell.col(1);
        cell.col(1) = d;
      }

      Lattice::t_Sites::iterator i_site = lattice->sites.begin(); 
      Lattice::t_Sites::iterator i_site_end = lattice->sites.end();
      for( types::t_unsigned n=0; i_site != i_site_end; ++i_site, ++n )
        i_site->site = n;

      if( not ( fill_structure( _structure ) ) ) return false;

      return true;
    }

  } // namespace Crystal

} // namespace LaDa
