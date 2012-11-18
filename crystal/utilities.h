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

#   ifdef LADA_STATIC
#     error LADA_STATIC already defined
#   endif
#   ifdef LADA_CRYSTAL_MODULE
#     define LADA_STATIC static
      //! Refolds a periodic vector into the unit cell.
      static math::rVector3d into_cell( math::rVector3d const &_vec, 
                                        math::rMatrix3d const &_cell, 
                                        math::rMatrix3d const &_inv);
      //! Refolds a periodic vector into the voronoi cell (eg first BZ or WignerSeitz).
      static math::rVector3d into_voronoi( math::rVector3d const &_vec, 
                                           math::rMatrix3d const &_cell, 
                                           math::rMatrix3d const &_inv);
      //! \brief Refolds a periodic vector into a cell centered around zero (in
      //!        fractional coordinates).
      //! \details Since the vector is refolded in fractional coordinates, it may
      //!          or may not be the vector with smallest norm. Use math::rVector3d
      //!          into_voronoi() to get the equivalent vector with smallest norm.
      static math::rVector3d zero_centered( math::rVector3d const &_vec, 
                                            math::rMatrix3d const &_cell, 
                                            math::rMatrix3d const &_inv);
#   else
#     define LADA_STATIC
      //! Refolds a periodic vector into the unit cell.
      inline math::rVector3d into_cell( math::rVector3d const &_vec, 
                                        math::rMatrix3d const &_cell, 
                                        math::rMatrix3d const &_inv)
        { return (*(math::rVector3d(*)(math::rVector3d const&, 
                                       math::rMatrix3d const&, 
                                       math::rMatrix3d const&))
                   api_capsule[24])(_vec, _cell, _inv); }
      //! Refolds a periodic vector into the voronoi cell (eg first BZ or WignerSeitz).
      inline math::rVector3d into_voronoi( math::rVector3d const &_vec, 
                                           math::rMatrix3d const &_cell, 
                                           math::rMatrix3d const &_inv)
        { return (*(math::rVector3d(*)(math::rVector3d const&, 
                                       math::rMatrix3d const&, 
                                       math::rMatrix3d const&))
                   api_capsule[25])(_vec, _cell, _inv); }
      //! Refolds a periodic vector into the voronoi cell (eg first BZ or WignerSeitz).
      //! \brief Refolds a periodic vector into a cell centered around zero (in
      //!        fractional coordinates).
      //! \details Since the vector is refolded in fractional coordinates, it may
      //!          or may not be the vector with smallest norm. Use math::rVector3d
      //!          into_voronoi() to get the equivalent vector with smallest norm.
      inline math::rVector3d zero_centered( math::rVector3d const &_vec, 
                                            math::rMatrix3d const &_cell, 
                                            math::rMatrix3d const &_inv)
        { return (*(math::rVector3d(*)(math::rVector3d const&, 
                                       math::rMatrix3d const&, 
                                       math::rMatrix3d const&))
                   api_capsule[26])(_vec, _cell, _inv); }
      //! Refolds a periodic vector into the voronoi cell (eg first BZ or WignerSeitz).
#   endif

    //! Refolds a periodic vector into the unit cell.
    LADA_STATIC inline math::rVector3d into_cell( math::rVector3d const &_vec, 
                                                  math::rMatrix3d const &_cell )
      { return into_cell(_vec, _cell, _cell.inverse()); }
    
   
    //! \brief Refolds a periodic vector into the voronoi cell (eg first BZ or WignerSeitz).
    //! \details May fail if the matrix is a weird parameterization of the
    //!          lattice. It is best to use a grubber(...) cell . 
    LADA_STATIC inline math::rVector3d into_voronoi( math::rVector3d const &_vec, 
                                                     math::rMatrix3d const &_cell )
      { return into_voronoi(_vec, _cell, _cell.inverse()); }

    //! \brief Refolds a periodic vector into a cell centered around zero (in
    //!        fractional coordinates).
    //! \details Since the vector is refolded in fractional coordinates, it may
    //!          or may not be the vector with smallest norm. Use math::rVector3d
    //!          into_voronoi() to get the equivalent vector with smallest norm.
    LADA_STATIC inline math::rVector3d zero_centered( math::rVector3d const &_vec, 
                                                      math::rMatrix3d const &_cell )
      { return zero_centered(_vec, _cell, _cell.inverse()); }

#   undef LADA_STATIC
  } // namespace Crystal
} // namespace LaDa
  
#endif
