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
    //! Creates divide and conquer box with periodic boundary condition.
    static PyObject* dnc_boxes_impl_( const Structure &_structure, 
                                      math::iVector3d const &_mesh, 
                                      types::t_real _overlap);

    namespace details
    {
      //! \brief Wrapper to python for periodic boundary divide and conquer.
      //! \see  LaDa::crystal::periodic_dnc()
      PyObject* pyperiodic_dnc(PyObject* _module, PyObject* _args, PyObject *_kwargs);
    }
  } // namespace crystal

} // namespace LaDa


#endif
