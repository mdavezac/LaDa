//
//  Version: $Id$
//
#include "va.h"

namespace Pescan
{ 

  types::t_real VirtualAtom :: potential_gradient( types::t_unsigned _pos )
  {
    t_Atom &atom = structure.atoms[ _pos ];
    Bands copy_bands = t_PescanBase::escan.bands; // last known result. copy em
    Bands result = copy_bands;

    //! flip atom
    atom->type = atom->type > 0 ? -1.0: 1.0; 
    t_PescanBase
  }

}
