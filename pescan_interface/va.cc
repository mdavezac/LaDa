//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "va.h"

namespace LaDa
{
  namespace Pescan
  { 
    VirtualAtom::t_Type VirtualAtom :: evaluate_one_gradient( types::t_unsigned _pos )
    {
      //! finds index of atom to change
      types::t_unsigned index = 0;
      t_Atoms :: const_iterator i_atom = structure.atoms.begin();
      t_Atoms :: const_iterator i_atom_end = structure.atoms.end();
      for(++_pos; i_atom != i_atom_end and _pos; ++i_atom, ++index )
        if( not ( i_atom->freeze & t_Atom::FREEZE_T ) ) --_pos;

      types::t_real result = 0;
      if( do_gradients & CHEMICAL_GRADIENT )
      {
        vff.chemical( index );
        result = apply_wfns() / t_Vff::deriv_amplitude;
      }
      if( do_gradients & STRESS_GRADIENT )
      {
        vff.stress( index );
        result += apply_wfns() / t_Vff::deriv_amplitude;
      }

      return result;
    }

  }
} // namespace LaDa
