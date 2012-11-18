#ifndef LADA_CRYSTAL_SPACEGROUP_H
#define LADA_CRYSTAL_SPACEGROUP_H

#include "LaDaConfig.h"

#include "structure/structure.h"

namespace LaDa
{
  namespace crystal 
  {
    namespace
    {
      //! \brief Finds and stores point group operations.
      //! \details Rotations are determined from G-vector triplets with the same
      //!          norm as the unit-cell vectors.
      //! \param[in] cell The cell for which to find the point group.
      //! \param[in] _tol acceptable tolerance when determining symmetries.
      //!             -1 implies that types::tolerance is used.
      //! \retval python list of affine symmetry operations for the given structure.
      //!         Each element is a 4x3 numpy array, with the first 3 rows
      //!         forming the rotation, and the last row is the translation.
      //!         The affine transform is applied as rotation * vector + translation.
      //!         `cell_invariants` always returns isometries (translation is zero).
      //! \see Taken from Enum code, PRB 77, 224115 (2008).
      PyObject* cell_invariants(math::rMatrix3d const &_cell, types::t_real _tolerance = -1e0)
#     ifdef LADA_CRYSTAL_MODULE
        ;
#     else
          { return (*(PyObject*(*)(math::rMatrix3d const &, types::t_real))
                    api_capsule[14])(_cell, _tolerance); }
#     endif
      //! \brief Finds and stores space group operations.
      //! \param[in] _structure The structure for which to find the space group.
      //! \param[in] _tol acceptable tolerance when determining symmetries.
      //!             -1 implies that types::tolerance is used.
      //! \retval spacegroup python list of symmetry operations for the given structure.
      //!         Each element is a 4x3 numpy array, with the first 3 rows
      //!         forming the rotation, and the last row is the translation.
      //!         The affine transform is applied as rotation * vector + translation.
      //! \warning Works for primitive lattices only.
      //! \see Taken from Enum code, PRB 77, 224115 (2008).
      PyObject* space_group(Structure const &_lattice, types::t_real _tolerance = -1e0)
#     ifdef LADA_CRYSTAL_MODULE
        ;
#     else
          { return (*(PyObject*(*)(Structure const &, types::t_real))
                    api_capsule[15])(_lattice, _tolerance); }
#     endif
    } // anonymous namespace
  }
}

#endif
