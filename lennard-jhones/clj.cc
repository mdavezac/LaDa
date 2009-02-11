//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "clj.h"

namespace LaDa
{
  //! Wrapper around Xiuwen's Lennard-Jhones+Coulomb functional.
  namespace CLJ
  {
    namespace details 
    {
      inline void anint( atat::rVector3d& d )
      {
        d[0] = std::floor( d[0] );
        d[1] = std::floor( d[1] );
        d[2] = std::floor( d[2] );
      }
    }

    Clj :: t_Return Clj :: all_( const t_Arg& _in, t_Arg &_out )
    {
      if( _in.atoms.size() != _out.atoms.size() ) 
        _out.atoms.resize( _in.atoms.size() );
      const atat::rMatrix3d inv_cell( !_in.cell );

      typedef TStructure ::  t_Atoms :: const_iterator t_cit;
      TStructure :: t_Atoms :: const_iterator i_atom_begin = _in.atoms.begin();
      TStructure :: t_Atoms :: const_iterator i_atom_end = _in.atoms.end();
      TStructure :: t_Atoms :: iterator i_force = _out.atoms.begin();
      for( t_cit i_atom1( i_atom_begin ); i_atom1 != i_atom_end; ++i_atom1 )
      {
        __DOASSERT( species_.end() != species_.find( i_atom1->type ),
                    "Specie " + i_atom1->type + " does not exist.\n" )

        const Specie &specie1( species_[ i_atom1->type ] );
        const atat::rVector3d fractional1( inv_cell * i_atom1->type );

        for( t_cit i_atom2( i_atom1 + 1 ); i_atom2 != i_atom_end; ++i_atom2 )
        {
          __DOASSERT( species_.end() != species_.find( i_atom2->type ),
                      "Specie " + i_atom2->type + " does not exist.\n" )
          const Specie &specie2( species_[ i_atom2->type ] );
          const Specie::t_Radius radius( specie1.radius + specie2.radius );
          const atat::rVector3d fractional2( inv_cell * i_atom1->type );
          const Specie :: t_Radius rcut( rcut_lj * radius * rcut_lj * radius );

          for( size_t i(-lj_mesh[0]); i < lj_mesh[0]; ++i )
            for( size_t j(-lj_mesh[1]); j < lj_mesh[1]; ++j )
              for( size_t k(-lj_mesh[2]); k < lj_mesh[2]; ++k )
              {
                // computes distance.
                atat::rVector3d d( fractional1 - fractional2 );
                details::anint( d );
                d += atat::rVector3d( i, j, k );
                d = _in.cell * d;
                if( atat::norm2( d) > rcut ) continue;

                single_( *i_force, _out.cell, d, rcut );
              }
        }
      }
    }

  } // namespace CLJ
} // namespace LaDa
