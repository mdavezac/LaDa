//
//  Version: $Id: vff.h 895 2008-12-22 02:04:18Z davezac $
//
#ifndef _MODELS_EWALD_FUNCTIONAL_H_
#define _MODELS_EWALD_FUNCTIONAL_H_

#include "PyladaConfig.h"

#include <Python.h>

namespace Pylada
{
  namespace pcm
  {
    //! Python C-interface to the fortran ewald function.
    PyObject* ewald(PyObject *_module, PyObject* _args, PyObject* _kwargs);
  } // namespace models.
} // namespace Pylada

#endif // _VFF_FUNCTIONAL_H_
