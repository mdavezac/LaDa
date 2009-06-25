//
//  Version: $Id$
//

#ifndef _LADA_CRYSTAL_SMITH_H_
#define _LADA_CRYSTAL_SMITH_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/tuple/tuple.hpp>

#include <atat/vectmac.h>
#include "structure.h"

namespace LaDa 
{
  namespace Crystal 
  {
    //! The type holding the smith transform.
    typedef boost::tuples::tuple<atat::rMatrix3d, atat::iVector3d> t_SmithTransform;

    //! Returns smith transform.
    t_SmithTransform get_smith_transform( atat::rMatrix3d const &_lat_cell,
                                          atat::rMatrix3d const &_str_cell );
    //! Returns smith transform.
    template< class T_TYPE >
      t_SmithTransform get_smith_transform( Crystal::TStructure<T_TYPE> const &_structure )
      {
        __DOASSERT( _structure.lattice == NULL, "Lattice not set in structure.\n" );
        return get_smith_transform( _structure.lattice->cell, _structure.cell ); 
      }

    //! Computes smith indices of position \a _pos.
    atat::iVector3d get_smith_index( t_SmithTransform const &_transformation,
                                     atat::rVector3d  const &_pos );
  } // namespace Crystal

} // namespace LaDa

#endif
