//
//  Version: $Id$
//

#ifndef _LADA_CRYSTAL_FRACTIONAL_CARTESIAN_H_
#define _LADA_CRYSTAL_FRACTIONAL_CARTESIAN_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "structure.h"

namespace LaDa 
{
  namespace Crystal {

    //! Converts a structure to fractional coordinates.
    template<class T_TYPE> 
      void to_fractional( TStructure<T_TYPE> &_structure );
    //! Converts a structure to cartesian coordinates.
    template<class T_TYPE> 
      void to_cartesian( TStructure<T_TYPE> &_structure );

    template<class T_TYPE> 
      void to_fractional( TStructure<T_TYPE> &_structure )
      {
        const atat::rMatrix3d inv( !_structure.cell );
        typename TStructure<T_TYPE>::t_Atoms::iterator i_atom = _structure.atoms.begin();
        typename TStructure<T_TYPE>::t_Atoms::iterator i_atom_end = _structure.atoms.end();
        for(; i_atom != i_atom_end; ++i_atom )
          i_atom->pos = inv * i_atom->pos;
      }
    template<class T_TYPE> 
      void to_cartesian( TStructure<T_TYPE> &_structure )
      {
        const atat::rMatrix3d &cell( _structure.cell );
        typename TStructure<T_TYPE>::t_Atoms::iterator i_atom = _structure.atoms.begin();
        typename TStructure<T_TYPE>::t_Atoms::iterator i_atom_end = _structure.atoms.end();
        for(; i_atom != i_atom_end; ++i_atom )
          i_atom->pos = cell * i_atom->pos;
      }
  } // namespace Crystal

} // namespace LaDa

#endif
