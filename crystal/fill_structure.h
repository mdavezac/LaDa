#ifndef _LADA_CRYSTAL_FILL_STRUCTURE_H_
#define _LADA_CRYSTAL_FILL_STRUCTURE_H_

#include "LaDaConfig.h"

#include "structure.h"

namespace LaDa 
{
  namespace Crystal 
  {
    //! \brief fills in \a atoms member of an Crystal::Structure instance from the
    //!        cell-shape and the lattice.
    bool fill_structure( Crystal::Structure &_str );
    //! \brief fills in \a atoms member of an Crystal::Structure instance from the
    //!        cell-shape and the lattice.
    bool fill_structure( Crystal::TStructure<std::string> &_str );
  } // namespace Crystal

} // namespace LaDa

#endif
