//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <algorithm>
#include <functional>

#include <opt/ndim_iterator.h>
#include <opt/atat.h>

#include "vff.h"
  
namespace LaDa
{
  namespace Vff
  { 
    bool Vff :: initialize_centers()
    {
      centers.clear();
      
      // Creates a list of centers
      t_Atoms :: iterator i_atom = structure.atoms.begin();
      t_Atoms :: iterator i_atom_end = structure.atoms.end();
      for(types::t_unsigned index=0; i_atom != i_atom_end; ++i_atom, ++index )
        centers.push_back( AtomicCenter( structure, *i_atom, index ) );

      // finds first neighbors on ideal lattice.
      typedef std::vector< std::vector< atat::rVector3d > > t_FirstNeighbors;
      t_FirstNeighbors fn;
      first_neighbors_( fn );


      // Creates a list of closest neighbors
      std::vector< atat::rVector3d > neighbors;
      typedef Crystal::Lattice::t_Sites :: iterator t_it;
      t_it i_site_begin = structure.lattice->sites.begin();
      t_it i_site, i_site2;
      t_it i_site_end = structure.lattice->sites.end();
      
      for(i_site = i_site_begin; i_site != i_site_end; ++i_site )
      {
        for(i_site2 = i_site_begin; i_site2 != i_site_end; ++i_site2 )
        {
          opt::Ndim_Iterator< types::t_int, std::less_equal<types::t_int> > period;
          period.add(-1,1);
          period.add(-1,1);
          period.add(-1,1);
          do // goes over periodic image
          {
            // constructs perdiodic image of atom *i_bond
            atat::rVector3d frac_image, image;
            frac_image[0] =  (types::t_real) period.access(0);
            frac_image[1] =  (types::t_real) period.access(1);
            frac_image[2] =  (types::t_real) period.access(2);
            image = i_site2->pos + structure.lattice->cell * frac_image;
            if( atat::norm2( image - i_site->pos ) > types::tolerance )
              neighbors.push_back( image - i_site->pos );
          } while ( ++period ); 
        }
      }

      // Sorts the neighbors according to distance from origin
      std::sort( neighbors.begin(), neighbors.end(), atat::norm_compare() );
      // And reduces to first neighbors only
      neighbors.resize(4*structure.lattice->sites.size());

      t_Centers :: iterator i_begin = centers.begin();
      t_Centers :: iterator i_end = centers.end();
      t_Centers :: iterator i_center, i_bond;
      atat::rVector3d frac_image;
      atat::rVector3d cut;
      cut = neighbors.front();
      types::t_real cutoff = types::t_real(0.25) * atat::norm2( neighbors.front() );
      for( i_center = i_begin; i_center != i_end; ++i_center )
      {
        for( i_bond = i_begin; i_bond != i_end; ++i_bond)
        {
          if( i_bond == i_center ) continue;
          
          std::vector<atat::rVector3d> :: const_iterator i_neigh = neighbors.begin();
          std::vector<atat::rVector3d> :: const_iterator i_neigh_end = neighbors.end();
          for(; i_neigh != i_neigh_end; ++i_neigh )
          {
            const atat::rVector3d image
            ( 
              i_center->origin->pos - *i_neigh - i_bond->origin->pos 
            );
            const atat::rVector3d frac_image( (!structure.cell) * image );
            const atat::rVector3d frac_centered
            ( 
              frac_image[0] - rint( frac_image[0] ),
              frac_image[1] - rint( frac_image[1] ),
              frac_image[2] - rint( frac_image[2] )
            );
            const atat::rVector3d cut( structure.cell * frac_centered );

            if( atat::norm2( cut ) > cutoff ) continue;
            
            i_center->bonds.push_back( t_Center ::__make__iterator__(i_bond) );
            const atat::rVector3d trans
            (
              rint( frac_image[0] ),
              rint( frac_image[1] ),
              rint( frac_image[2] ) 
            );
            i_center->translations.push_back( trans );
            i_center->do_translates.push_back
            ( 
              atat::norm2(trans) > atat::zero_tolerance 
            );

            if( i_center->bonds.size() == 4 ) break;
          } // loop over neighbors
          if( i_center->bonds.size() == 4 ) break;
        } // loop over bonds
      } // loop over centers 

      __DODEBUGCODE( check_tree(); )
      return true;
    } // Vff :: construct_bond_list

  } // namespace vff
} // namespace LaDa
