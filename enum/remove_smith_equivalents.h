#ifndef LADA_ENUM_REMOVE_TRANSLATIONAL_EQUIVALENTS_H_
#define LADA_ENUM_REMOVE_TRANSLATIONAL_EQUIVALENTS_H_

#include "LaDaConfig.h"

#include "numeric_type.h"

namespace LaDa
{
  //! \cond
  namespace Crystal { class Lattice; } 
  //! \endcond

  namespace enumeration
  {
    //! \cond
    class SmithGroup;
    //! \endcond

    //! Turns off translational equivalents in a SmithGroup.
    void remove_smith_equivalents( Database &_database,
                                   SmithGroup const &_sg, 
                                   Crystal::Lattice const &_lattice );
  }
}

#endif
