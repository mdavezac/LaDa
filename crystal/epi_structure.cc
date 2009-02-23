//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/debug.h>
#include <opt/types.h>

#include "lattice.h"
#include "structure.h"
#include "epi_structure.h"
#include "fill_structure.h"



namespace LaDa
{

  namespace Crystal 
  {
    bool create_epitaxial_structure( Structure& _structure,
                                     atat::rVector3d &_direction,
                                     atat::iVector3d &_extent )
    {

      // Then constructs unit cell
      atat::rMatrix3d &cell( _structure.cell ); 
      atat::rMatrix3d diagonal;
      diagonal.set_diagonal( (types::t_real) _extent(0),
                             (types::t_real) _extent(1),
                             (types::t_real) _extent(2)  );
      Lattice *lattice( _structure.lattice );
      cell = lattice->cell;
      cell.set_column(0, lattice->cell * _direction ); 
      cell = cell * diagonal;

      // Checks that cell is not singular
      if ( Fuzzy::is_zero( std::abs(det(cell)) ) )
        cell.set_column(1, lattice->cell.get_column( 0 ) );
      if ( Fuzzy::is_zero( std::abs(det(cell)) ) )
        cell.set_column(2, lattice->cell.get_column( 1 ) );
      if ( Fuzzy::is_zero( std::abs(det(cell)) ) )
      {
        std::cerr << "Could not construct unit-cell\n" << cell << std::endl;
        return false;
      }
      if( Fuzzy::is_zero( det(cell) ) )
      {
        atat::rVector3d swap = cell. get_column(1);
        cell.set_column(1, cell.get_column( 2 ));
        cell.set_column(2, swap);
      }

      // Makes sure the triad is direct
      if ( det(cell) < 0 )
      {
        atat::rVector3d d = cell.get_column(2);
        cell.set_column(2, cell.get_column(1) );
        cell.set_column(1, d);
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
