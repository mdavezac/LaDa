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

    //! Computes linear smith index of position \a _pos.
    inline size_t get_linear_smith_index( t_SmithTransform const &_transformation,
                                          atat::rVector3d  const &_pos )
    {
      atat::iVector3d const indices = get_smith_index( _transformation, _pos );
      atat::iVector3d const &smith = boost::tuples::get<1>(_transformation);
      return indices(2) + smith(2) * ( indices(1) + indices(0) * smith(1) );
    }


    inline bool equivalent_mod_cell( t_SmithTransform const &_transformation, 
                                     atat::rVector3d const &_a, atat::rVector3d const &_b )
      { return get_smith_index(_transformation, _a - _b) == atat::iVector3d(0,0,0); }

    //! Computes map of smith indices.
    void get_smith_map( Crystal::Structure const &_str, std::vector< std::vector<size_t> > &_out);
  } // namespace Crystal

} // namespace LaDa

#endif
