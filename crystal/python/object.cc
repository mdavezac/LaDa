#include <LaDaConfig.h>

#include "object.h"

namespace LaDa 
{
  namespace python
  {
    namespace 
    {
      void object_reset(Object& _self, PyObject *_in = NULL)
      {
        PyObject * const dummy(_self.object_);
        _self.object_ = _in;
        Py_XINCREF(_self.object_);
        Py_XDECREF(dummy);
      }
    }
  }
}
