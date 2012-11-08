#ifndef LADA_CRYSTAL_UTILITIES_H
#define LADA_CRYSTAL_UTILITIES_H

#include "LaDaConfig.h"


#include <misc/types.h>
#include <math/misc.h>
#include "atom/atom.h"

namespace LaDa
{
  namespace crystal 
  {

    //! Refolds a periodic vector into the unit cell.
    math::rVector3d into_cell( math::rVector3d const &_vec, 
                               math::rMatrix3d const &_cell, 
                               math::rMatrix3d const &_inv);
    //! Refolds a periodic vector into the unit cell.
    inline math::rVector3d into_cell( math::rVector3d const &_vec, 
                               math::rMatrix3d const &_cell )
      { return into_cell(_vec, _cell, _cell.inverse()); }
    
   
    //! Refolds a periodic vector into the voronoi cell (eg first BZ or WignerSeitz).
    math::rVector3d into_voronoi( math::rVector3d const &_vec, 
                                  math::rMatrix3d const &_cell, 
                                  math::rMatrix3d const &_inv);
    //! \brief Refolds a periodic vector into the voronoi cell (eg first BZ or WignerSeitz).
    //! \details May fail if the matrix is a weird parameterization of the
    //!          lattice. It is best to use a grubber(...) cell . 
    inline math::rVector3d into_voronoi( math::rVector3d const &_vec, 
                                  math::rMatrix3d const &_cell )
      { return into_voronoi(_vec, _cell, _cell.inverse()); }

    //! \brief Refolds a periodic vector into a cell centered around zero (in
    //!        fractional coordinates).
    //! \details Since the vector is refolded in fractional coordinates, it may
    //!          or may not be the vector with smallest norm. Use math::rVector3d
    //!          into_voronoi() to get the equivalent vector with smallest norm.
    math::rVector3d zero_centered( math::rVector3d const &_vec, 
                                   math::rMatrix3d const &_cell, 
                                   math::rMatrix3d const &_inv);
    //! \brief Refolds a periodic vector into a cell centered around zero (in
    //!        fractional coordinates).
    //! \details Since the vector is refolded in fractional coordinates, it may
    //!          or may not be the vector with smallest norm. Use math::rVector3d
    //!          into_voronoi() to get the equivalent vector with smallest norm.
    inline math::rVector3d zero_centered( math::rVector3d const &_vec, 
                                   math::rMatrix3d const &_cell )
      { return zero_centered(_vec, _cell, _cell.inverse()); }

  } // namespace Crystal
} // namespace LaDa
  
#endif
