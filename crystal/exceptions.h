#ifndef LADA_CRYSTAL_EXCEPTIONS_H
#define LADA_CRYSTAL_EXCEPTIONS_H

#include "LaDaConfig.h"

#include <iostream>

#include <math/exceptions.h>


namespace LaDa
{
  namespace error
  {
    //! Root of execptions thrown in crystal namespace.
    struct crystal: virtual root {};
    //! Thrown when a function requires the structure have only one specie per atom.
    struct too_many_species: virtual crystal, virtual input {};
    //! Thrown when structure does not contain atomic sites.
    struct empty_structure: virtual crystal, virtual input {};

    //! Thrown when the number of atoms in a supercell does not match the lattice.
    struct incommensurate_number_of_sites: virtual crystal, virtual unideal_lattice, virtual input {};
    //! Thrown when the site index has not been initialized.
    struct incorrect_site_index: virtual crystal, virtual unideal_lattice {};
    //! Thrown when a lattice is not primitive.
    struct not_primitive : virtual crystal, virtual input {};
  }
}

#endif
