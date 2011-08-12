#ifndef LADA_MATH_EXCEPTIONS_H_
#define LADA_MATH_EXCEPTIONS_H_

#include "LaDaConfig.h"

#include <root_exceptions.h>

namespace LaDa
{
  namespace error
  {
    //! Root of input errors.
    struct array_of_different_sizes: virtual internal {};
    //! \brief Thrown when an structure is not the expected supercell of a lattice.
    //! \details This error may occur in methods which expects structures to be
    //!          supercells of a lattice, without allowance for relaxation. 
    struct unideal_lattice: virtual root {};
    //! Thrown when an atomic position is unexpectedly off-lattice.
    struct off_lattice_position: virtual unideal_lattice {};
    //! Thrown when a structure is not a supercell of an ideal lattice.
    struct not_a_supercell: virtual unideal_lattice, virtual input {};
  }
}

#endif
