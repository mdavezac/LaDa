#ifndef LADA_CRYSTAL_EXCEPTIONS_H
#define LADA_CRYSTAL_EXCEPTIONS_H

#include "LaDaConfig.h"

#include <iostream>

#include <math/exceptions.h>


namespace LaDa
{
  namespace error
  {
    //! Thrown when a function requires the structure have only one specie per atom.
    struct too_many_species: virtual input {};
    //! Thrown when structure does not contain atomic sites.
    struct empty_structure: virtual input {};

    //! Thrown when the number of atoms in a supercell does not match the lattice.
    struct incommensurate_number_of_sites: virtual unideal_lattice, virtual input {};
    //! Thrown when the site index has not been initialized.
    struct incorrect_site_index: virtual unideal_lattice {};
  }
}

#endif
