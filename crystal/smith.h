//
//  Version: $Id$
//

#ifndef _LADA_CRYSTAL_SMITH_H_
#define _LADA_CRYSTAL_SMITH_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/tuple/tuple.hpp>

#include "structure.h"

namespace LaDa 
{
  namespace Crystal 
  {
    //! The type holding the smith transform.
    typedef boost::tuples::tuple<math::rMatrix3d, math::iVector3d> t_SmithTransform;

    //! Returns smith transform.
    t_SmithTransform get_smith_transform( math::rMatrix3d const &_lat_cell,
                                          math::rMatrix3d const &_str_cell );
    //! Returns smith transform.
    template< class T_TYPE >
      t_SmithTransform get_smith_transform( Crystal::TStructure<T_TYPE> const &_structure )
      {
        __DOASSERT( _structure.lattice == NULL, "Lattice not set in structure.\n" );
        return get_smith_transform( _structure.lattice->cell, _structure.cell ); 
      }

    //! Computes smith indices of position \a _pos.
    math::iVector3d get_smith_index( t_SmithTransform const &_transformation,
                                     math::rVector3d  const &_pos );

    //! Computes linear smith index from non-linear smith index.
    inline size_t get_linear_smith_index( math::iVector3d const &_extent,
                                          math::iVector3d  const &_index )
      { return _index(2) + _extent(2) * ( _index(1) + _index(0) * _extent(1) ); }

    //! Computes linear smith index from non-linear smith index, including sites.
    inline size_t get_linear_smith_index( math::iVector3d const &_extent,
                                          size_t const &_site_index,
                                          math::iVector3d  const &_index )
      { return _index(2)+_extent(2)*(_index(1)+_extent(1)*(_index(0)+_site_index*_extent(0))); }

    //! Computes linear smith index of position \a _pos.
    inline size_t get_linear_smith_index( t_SmithTransform const &_transformation,
                                          math::rVector3d  const &_pos )
    {
      return get_linear_smith_index
             (
               boost::tuples::get<1>(_transformation),
               get_smith_index( _transformation, _pos )
             );
    }
    //! Computes linear smith index of position \a _pos.
    inline size_t get_linear_smith_index( t_SmithTransform const &_transformation,
                                          size_t const &_site_index,
                                          math::rVector3d  const &_pos )
    {
      return get_linear_smith_index
             (
               boost::tuples::get<1>(_transformation),
               _site_index,
               get_smith_index( _transformation, _pos )
             );
    }
    //! Computes linear smith index of position \a _pos.
    template<class T_TYPE>
      size_t get_linear_smith_index( t_SmithTransform const &_transformation,
                                     Crystal::Atom_Type<T_TYPE> &_atom ) 
      {
        LADA_ASSERT( _atom.site != -1, "Site index not set.\n");
        return get_linear_smith_index
               (
                 boost::tuples::get<1>(_transformation),
                 _atom.site,
                 get_smith_index( _transformation, _atom.pos )
               );
      }

    inline bool are_periodic_images( t_SmithTransform const &_transformation, 
                                     math::rVector3d const &_a, math::rVector3d const &_b )
      { return get_smith_index(_transformation, _a - _b) == math::iVector3d(0,0,0); }

    //! Computes map of smith indices.
    void get_smith_map( Crystal::Structure const &_str, std::vector< std::vector<size_t> > &_out);
  } // namespace Crystal

} // namespace LaDa

#endif
