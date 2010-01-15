//
//  Version: $Id$
//

#ifndef _LADA_EPI_STRUCTURE_H_
#define _LADA_EPI_STRUCTURE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <Eigen/Core>

namespace LaDa 
{
  namespace Crystal
  {
    //! \cond
    class Structure;
    //! \endcond

    //! \brief Creates an epitaxial structure.
    //! \param[out] _structure the structure on output. The atoms are ordered
    //!                        with respect to the growth direction.
    //! \param[in] _direction the growth direction.
    //! \param[in] _extent the size in the direction growth ( x ), the other two
    //!                    directions. These are taken from the lattice
    //!                    unit-cell such that the determinant of the structure
    //!                    cell is strictly positive. The first vector in the
    //!                    unit-cell is the growth direction.
    bool create_epitaxial_structure( Structure& _structure,
                                     Eigen::Vector3d &_direction,
                                     Eigen::Vector3i &_extent );

  } // namespace Crystal

} // namespace LaDa

#endif
