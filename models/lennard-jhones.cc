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
      if( _in.atoms.size() != _out.atoms.size() ) 
        _out.atoms.resize( _in.atoms.size() );
      const atat::rMatrix3d inv_cell( !_in.cell );
      types::t_real energy(0);

      typedef TStructure ::  t_Atoms :: const_iterator t_cit;
      TStructure :: t_Atoms :: const_iterator i_atom_begin = _in.atoms.begin();
      TStructure :: t_Atoms :: const_iterator i_atom_end = _in.atoms.end();
      TStructure :: t_Atoms :: iterator i_force1 = _out.atoms.begin();
      for( t_cit i_atom1( i_atom_begin ); i_atom1 != i_atom_end; ++i_atom1, ++i_force1 )
      {
        __DOASSERT( species_.end() != species_.find( i_atom1->type ),
                    "Specie " + i_atom1->type + " does not exist.\n" )

        const Specie &specie1( species_[ i_atom1->type ] );

        TStructure :: t_Atoms :: iterator i_force2 = i_force1 + 1;
        for( t_cit i_atom2( i_atom1 + 1 ); i_atom2 != i_atom_end; ++i_atom2, ++i_force2 )
        {
          __DOASSERT( species_.end() != species_.find( i_atom2->type ),
                      "Specie " + i_atom2->type + " does not exist.\n" )
          const Specie &specie2( species_[ i_atom2->type ] );
          const Specie::t_Radius sigma( specie1.radius + specie2.radius );
          const Specie :: t_Radius sigma_squared( sigma * sigma );
          const Specie :: t_Radius rcut_squared( rcut * rcut * sigma_squared );

          const atat::rVector3d dfractional
                                ( 
                                  std::floor( i_atom1->pos[0] - i_atom2->pos[0] ),
                                  std::floor( i_atom1->pos[1] - i_atom2->pos[1] ),
                                  std::floor( i_atom1->pos[2] - i_atom2->pos[2] )
                                );
          for( size_t i(-mesh[0]); i < mesh[0]; ++i )
            for( size_t j(-mesh[1]); j < mesh[1]; ++j )
              for( size_t k(-mesh[2]); k < mesh[2]; ++k )
              {
                // computes distance.
                const atat::rVector3d distance( _in.cell * ( dfractional + atat::rVector3d(i,j,k) ) );
                const types::t_real normd( atat::norm2(d) );
                if( normd > rcut_squared ) continue;

                // energy -= 4.0 * scale * lj_pot( sigma_squared / rcut_squared )
                { // compute correction energy
                  const types::t_real sqared( sigma_squared / rcut_squared )
                  const types::t_real sixth( squared * squared *squared )
                  const types::t_real twelveth( squared * squared *squared )
                  energy -= 4.e0 * bond_strength * ( twelveth - sixth );
                }

                // energy += 4.0 * scale * lj_pot( sigma_squared / rcut_squared )
                // force += 24.0 * scale * ( 2*u^14 - u^8 ) * d
                // stress += force_alpha distance_beta
                { // compute correction energy
                  const types::t_real inv_normd( 1e0 / normd )
                  const types::t_real sqared( sigma_squared * inv_normd )
                  const types::t_real sixth( squared * squared *squared )
                  const types::t_real twelveth( squared * squared *squared )
                  energy -= 4.e0 * bond_strength * ( twelveth - sixth );

                  const types::t_real ffactor( 24.e0 * scale * ( 2.e0 * twelveth - sixth ) * inv_normd );
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
