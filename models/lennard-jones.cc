//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "lennard-jones.h"

namespace LaDa
{
  namespace Models
  {
    LennardJones :: t_Return LennardJones :: energy( const t_Arg& _in, t_Arg &_out )
    {
      __DOASSERT( _in.atoms.size() != _out.atoms.size(), "Incoherent structure size.\n" )
      types::t_real energy(0);

      typedef t_Arg ::  t_Atoms :: const_iterator t_cit;
      t_cit i_atom_begin = _in.atoms.begin();
      t_cit i_atom_end = _in.atoms.end();
      t_Arg :: t_Atoms :: iterator i_force1 = _out.atoms.begin();
      for( t_cit i_atom1( i_atom_begin ); i_atom1 != i_atom_end; ++i_atom1, ++i_force1 )
      {
        t_Arg :: t_Atoms :: iterator i_force2 = i_force1 + 1;
        for( t_cit i_atom2( i_atom1 + 1 ); i_atom2 != i_atom_end; ++i_atom2, ++i_force2 )
        {
          __DOASSERT( species_.end() != species_.find( i_atom1->type + i_atom2->type ),
                      "Bond " + i_atom1->type + "-" + i_atom2->type + " does not exist.\n" )
          const Bond &bond( species_[ i_atom2->type ] );
          const types::t_real rcut_squared( rcut_ * rcut_ );

          const atat::rVector3d dfractional
                                ( 
                                  std::floor( i_atom1->pos[0] - i_atom2->pos[0] ),
                                  std::floor( i_atom1->pos[1] - i_atom2->pos[1] ),
                                  std::floor( i_atom1->pos[2] - i_atom2->pos[2] )
                                );
          for( size_t i(-mesh_[0]); i < mesh_[0]; ++i )
            for( size_t j(-mesh_[1]); j < mesh_[1]; ++j )
              for( size_t k(-mesh_[2]); k < mesh_[2]; ++k )
              {
                // computes distance.
                const atat::rVector3d distance
                (
                  _in.cell * ( dfractional + atat::rVector3d(i,j,k) ) 
                );
                const types::t_real normd( atat::norm2(d) );
                if( normd > rcut_squared ) continue;

                // energy -= 4.0 * scale * lj_pot( sigma_squared / rcut_squared )
                { // compute cutoff correction energy
                  const types::t_real sqared( 1e0 / rcut_squared )
                  const types::t_real sixth( squared * squared * squared )
                  const types::t_real twelveth( sixth * sixth );
                  energy -=  bond.hard_sphere * twelveth - bond.vand_der_walls * sixth;
                }

                { // van_der_walls energy, force, stress.
                  const types::t_real squared( 1e0 / normd )
                  const types::t_real sixth( squared * squared * squared )
                  const types::t_real twelveth( sixth * sixth )
                  energy -= bond.hard_sphere * twelveth - bond.van_der_walls * sixth;

                  const types::t_real ffactor
                  ( 
                    squared * ( 2.e0 * bond.hard_sphere * twelveth - bond.van_der_walls * sixth )
                  );
                  const atat::rVector3d force( ffactor * distance );
                  i_force1->pos -= force;
                  i_force2->pos += force;

                  for( size_t i(0); i < 3; ++i )
                    for( size_t j(0); j < 3; ++j )
                      _out.cell(i,j) += -force[i] * distance[j];
                }

              } // loop over periodic images.
        } // loop over atom 2
      } // loop over atom 1
    }

  } // namespace CLJ
} // namespace LaDa
