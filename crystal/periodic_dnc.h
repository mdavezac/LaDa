#ifndef LADA_CRYSTAL_PERIODIC_DNC_H
#define LADA_CRYSTAL_PERIODIC_DNC_H

#include "LaDaConfig.h"

#include <vector>

#include <math/eigen.h>

#include "structure/structure.h"


namespace LaDa 
{
  namespace crystal 
  {
    namespace
    {
      //! Creates divide and conquer box with periodic boundary condition.
      static PyObject* dnc_boxes( const Structure &_structure, 
                                  math::iVector3d const &_mesh, 
                                  types::t_real _overlap)
#     ifdef LADA_CRYSTAL_MODULE
          ;
#     else
          { return (*(PyObject*(*)(Structure const&, math::iVector3d const&,
                                   types::t_real))
                    api_capsule[19])(_structure, _mesh, _overlap); }
#     endif
    } // anonymous namespace
  } // namespace crystal
} // namespace LaDa


#endif
