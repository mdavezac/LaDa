#ifndef LADA_CRYSTAL_COORDINATION_SHELLS_H
#define LADA_CRYSTAL_COORDINATION_SHELLS_H

#include "LaDaConfig.h"

#include <boost/lambda/lambda.hpp>

#include <opt/types.h>
#include "structure/structure.h"

namespace LaDa
{
  namespace crystal
  {
    //! \brief Creates list of coordination shells up to given order.
    //! \returns A list of lists of tuples. The outer list is over coordination shells.
    //!          The inner list references the atoms in a shell.
    //!          Each innermost tuple contains a reference to the atom in question,
    //!          a translation vector to its periodic image inside the relevant shell, 
    //!          and the distance from the center to the relevant periodic image.
    //! \param[in] _structure : Structure for which to determine coordination shells.
    //! \param[in] _nshells : Number of shells to compute.
    //! \param[in] _center : Center of the coordination shells.
    //! \param[in] _tolerance : criteria to judge when a shell ends.
    //! \param[in] _natoms : Total number of neighbors to consider. Defaults to fcc + some security.
    PyObject* coordination_shells( crystal::Structure const &_structure, Py_ssize_t _nshells, 
                                   math::rVector3d const &_center,
                                   types::t_real _tolerance=types::tolerance,
                                   Py_ssize_t _natoms = 0 );

  } // end of crystal namespace.
} // namespace LaDa

#endif
