//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "ewald.h"

extern "C" void FC_FUNC( ewaldf, EWALDF )
                (
                  const double *const, // verbosity
                  const double *const, // Energy
                  const double *const, // forces (reduced)
                  const double *const, // forces (cartesian)
                  const double *const, // stress
                  const double *const, // number of atoms
                  const double *const, // reduced atomic coordinates.
                  const double *const, // atomic charges
                  const double *const, // cell vectors
                  const double *const  // dimension of arrays.
                );
namespace LaDa
{
  namespace Models
  {
    Ewald :: t_Return Ewald :: energy( const t_Arg& _in, t_Arg &_out )
    {
      __DOASSERT( _in.atoms.size() != _out.atoms.size(), "Incoherent structure size.\n" )
      types::t_real energy(0);

      const size_t natom( _in.atoms.size() );
      double charges[ natom ];
      double positions[ natom * 3 ], forces[ natoms * 3 ], cforces[ natoms * 3 ];
      double cell[ 9 ], stress[ 6 ];
      const double verbosity(0);
      double energy(0);
      const double n( natom );

      typedef  ::  t_Atoms :: const_iterator t_cit;
      t_cit i_atom = _in.atoms.begin();
      for( size_t i(0); i < natom; ++i, ++i_atom )
      {
        __ASSERT( charges_.find( i_atom->type ) == charges_.end(),
                  "Atomic charge does not exist.\n" )
        charges[i] = charges_[i_atom->type];
        positions[ i*3 ]     = i_atom->pos[0];
        positions[ i*3 + 1 ] = i_atom->pos[1];
        positions[ i*3 + 2 ] = i_atom->pos[2];
      } // loop over atoms

      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
          cell[ i + j*3 ] = _in.cell(i,j);

      FC_FUNC( ewaldf, EWALDF )
      (
        &verbosity,   // verbosity
        &energy,      // Energy
        forces,       // forces (reduced)
        cforces,      // forces (cartesian)
        stress,       // stress
        &n,           // number of atoms
        positions,    // reduced atomic coordinates.
        charges,      // atomic charges
        cell,         // cell vectors
        &n            // dimension of arrays.
      );
      // Copy (reduced) forces.
      t_Arg :: t_Atom :: iterator i_force = _out.atoms.begin();
      for( size_t i(0); i < natom; ++i, ++i_force )
      {
        i_force->pos[0] += forces[ i*3 ];
        i_force->pos[1] += forces[ i*3 + 1 ];
        i_force->pos[2] += forces[ i*3 + 2 ];
      } // loop over atoms
      // copy stress.
      _out.cell(0,0) += stress[0];
      _out.cell(1,1) += stress[1];
      _out.cell(2,2) += stress[2];
      _out.cell(0,1) += stress[3];
      _out.cell(1,0) += stress[3];
      _out.cell(1,2) += stress[4];
      _out.cell(2,1) += stress[4];
      _out.cell(0,2) += stress[5];
      _out.cell(2,0) += stress[5];
      
      return energy;
    }

  } // namespace CLJ
} // namespace LaDa
