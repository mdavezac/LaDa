#ifndef LADA_OPT_SMITH_NORMAL_FORM_H
#define LADA_OPT_SMITH_NORMAL_FORM_H

#include "LaDaConfig.h"

#include <boost/tuple/tuple.hpp>

#include "eigen.h"
#include "misc.h"

namespace LaDa 
{
  namespace math
  {
    //! The type holding the smith transform.
    typedef boost::tuples::tuple<math::rMatrix3d, math::iVector3d> t_SmithTransform;

    //! Computes smith normal form of a matrix \a  _S = \a _L \a _M \a _R.
    void smith_normal_form( math::iMatrix3d& _S, math::iMatrix3d & _L,
                            const math::iMatrix3d& _M, math::iMatrix3d &_R );

    //! Computes smith transform between supercell and unit cell.
    t_SmithTransform smith_transform( math::rMatrix3d const &_unitcell,
                                      math::rMatrix3d const &_supercell );

    //! Computes smith indices of position \a _pos.
    math::iVector3d smith_index( t_SmithTransform const &_transformation,
                                 math::rVector3d  const &_pos );

    //! \brief Computes linear smith index from non-linear smith index.
    //! \warning Only valid for lattices with a single site.
    inline size_t linear_smith_index( math::iVector3d const &_extent,
                                      math::iVector3d  const &_index )
      { return _index(2) + _extent(2) * ( _index(1) + _index(0) * _extent(1) ); }

    //! Computes linear smith index from non-linear smith index, including sites.
    inline size_t linear_smith_index( math::iVector3d const &_extent,
                                      size_t const &_site_index,
                                      math::iVector3d  const &_index )
      { return _index(2)+_extent(2)*(_index(1)+_extent(1)*(_index(0)+_site_index*_extent(0))); }

    //! Computes linear smith index of position \a _pos.
    inline size_t linear_smith_index( t_SmithTransform const &_transformation,
                                      math::rVector3d  const &_pos )
    {
      return linear_smith_index
             (
               boost::tuples::get<1>(_transformation),
               smith_index( _transformation, _pos )
             );
    }
    //! Computes linear smith index of position \a _pos.
    inline size_t linear_smith_index( t_SmithTransform const &_transformation,
                                      size_t const &_site_index,
                                      math::rVector3d  const &_pos )
    {
      return linear_smith_index
             (
               boost::tuples::get<1>(_transformation),
               _site_index,
               smith_index( _transformation, _pos )
             );
    }

    //! True if two vectors are periodic images of each other.
    inline bool are_periodic_images( t_SmithTransform const &_transformation, 
                                     math::rVector3d const &_a, math::rVector3d const &_b )
      { return eq(smith_index(_transformation, _a - _b), math::iVector3d(0,0,0)); }
  } // namespace opt

} // namespace LaDa

#endif
