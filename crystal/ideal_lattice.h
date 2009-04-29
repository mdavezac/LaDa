//
//  Version: $Id$
//

#ifndef _LADA_CRYSTAL_IDEAL_LATTICE_H_
#define _LADA_CRYSTAL_IDEAL_LATTICE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include <boost/tuple/tuple.hpp>

#include <opt/types.h>

#include "lattice.h"
#include "structure.h"


namespace LaDa 
{
  namespace Crystal 
  {
    //! \brief Retrieves the deformation to go from \a _structure to its
    //!        corresponding ideal on a lattice.
    //! The algorithm operates by searching for the \_nneighs first neighbors
    //! of the first atom, ! both on the ideal lattice and on the deformed
    //! structure. It then computes the deformation from a least square fit
    //! between  ideal and deformed positions: \f$\sum_{i,j}
    //! ||(1+\bar{\epsilon})R_i + t_i - R_j^{(0)}||^2 \f$, where \f$t_i\f$ is a
    //! translation.
    boost::tuples::tuple< atat::rMatrix3d, atat::rVector3d >
      retrieve_deformation( const Structure &_structure,
                            const size_t _nneigs = 12 );

    //! Reduces list of positions to first neighbors of first position.
    void find_first_neighbors( std::vector< atat::rVector3d > &_positions,
                               const atat::rMatrix3d &_cell,
                               const size_t _n );
  } // namespace Crystal

} // namespace LaDa

#endif
