#ifndef LADA_CRYSTAL_QUANTITY_H
#define LADA_CRYSTAL_QUANTITY_H

#include "LaDaConfig.h"
#include <Python.h>

#include <misc/types.h>


namespace LaDa
{
  namespace math
  {
    //! Checks whether this a quantity object.
    bool PyQuantity_Check(PyObject *_in);
    //! Returns a quantity object from a number and a unit in string.
    PyObject *PyQuantity_FromC(types::t_real const &_double, std::string const &_units);
    //! \brief Returns a quantity object from a number and a quantity object.
    //! \details if _number has itself a unit, then it should be convertible to
    //! _units.  However, in that case, _number is returned as is
    //! (Py_INCREF'ed).
    PyObject *PyQuantity_FromPy(PyObject *_number, PyObject *_units);
    //! \brief Returns as a number in specified units.
    //! \details If the number is 0, then one should check whether exception
    //! was thrown.
    types::t_real PyQuantity_Get(PyObject *_number, std::string const &_units);
  }
}
#endif
