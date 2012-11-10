#include "LaDaConfig.h"

#include <python/object.h>
#include <python/exceptions.h>
#include "quantity.h"


namespace LaDa
{
  namespace math
  {
    PyObject* UnitQuantityClass()
    {
      python::Object quant( PyImport_ImportModule("quantities") );
      if(not quant) return NULL;
      return PyObject_GetAttrString(quant.borrowed(), "UnitQuantity");
    }
    PyObject* QuantityClass()
    {
      python::Object quant( PyImport_ImportModule("quantities") );
      if(not quant) return NULL;
      return PyObject_GetAttrString(quant.borrowed(), "Quantity");
    }
    bool PyQuantity_Check(PyObject *_in)
    {
      if(not _in) return false;
      PyObject* quantity_class = QuantityClass();
      if(not quantity_class) 
      {
        PyErr_Clear();
        return false;
      }
      bool const resultA = PyObject_IsInstance(_in, quantity_class) == 1;
      Py_DECREF(quantity_class);
      if(resultA) return true;
      quantity_class = UnitQuantityClass();
      if(not quantity_class) 
      {
        PyErr_Clear();
        return false;
      }
      bool const resultB = PyObject_IsInstance(_in, quantity_class) == 1;
      Py_DECREF(quantity_class);
      return resultB;
    }

    PyObject *PyQuantity_FromC(types::t_real const &_double, std::string const &_units)
    {
      // creates global/local dictionary in order to run code.
      python::Object const globals = PyImport_ImportModule("quantities");
      if(not globals) return NULL;
      python::Object locals = PyDict_New();
      if(not locals) return NULL;
      python::Object number = PyFloat_FromDouble(_double);
      if(not number) return NULL;
      python::Object units = PyString_FromString(_units.c_str());
      if(not units) return NULL;
      if(PyDict_SetItemString(locals.borrowed(), "number", number.borrowed()) < 0)
        return NULL;
      if(PyDict_SetItemString(locals.borrowed(), "units", units.borrowed()) < 0)
        return NULL;
      python::Object result(PyRun_String( "quantity.Quantity(number, units)", Py_eval_input,
                                  PyModule_GetDict(globals.borrowed()), 
                                  locals.borrowed() ));
      if(not result) return NULL;
      return result.release();
    }

    bool PyQuantity_Convertible(PyObject *_a, PyObject *_b)
    {
      if(not PyQuantity_Check(_a)) return false;
      if(not PyQuantity_Check(_b)) return false;
      char rescale_str[] = "rescale";
      char s_str[] = "O";
      PyObject *result = PyObject_CallMethod(_a, rescale_str, s_str, _b);
      if(PyErr_Occurred())
      {
        PyErr_Clear();
        Py_XDECREF(result);
        return false;
      }
      Py_XDECREF(result);
      return true;
    }
    PyObject *PyQuantity_FromPy(PyObject *_number, PyObject *_units)
    {
      if(not PyQuantity_Check(_units))
      {
        LADA_PYERROR(TypeError, "Expected a quantities object.");
        return NULL;
      }
      if(PyQuantity_Check(_number))
      {
        if(not PyQuantity_Convertible(_number, _units))
        {
          LADA_PYERROR(TypeError, "Input quantities are not convertible.");
          return NULL;
        }
        Py_INCREF(_number);
        return NULL;
      }
      
      // creates global/local dictionary in order to run code.
      python::Object const globals = PyImport_ImportModule("quantities");
      if(not globals) return NULL;
      python::Object locals = PyDict_New();
      if(not locals) return NULL;
      if(PyDict_SetItemString(locals.borrowed(), "number", _number) < 0)
        return NULL;
      if(PyDict_SetItemString(locals.borrowed(), "unitdims", _units) < 0)
        return NULL;
      python::Object result(PyRun_String( "number * unitdims.units", Py_eval_input,
                                  PyModule_GetDict(globals.borrowed()), 
                                  locals.borrowed() ));
      if(not result) return NULL;
      return result.release();
    }


    types::t_real PyQuantity_Get(PyObject *_in, std::string const &_units)
    {
      if(not PyQuantity_Check(_in))
      {
        if(PyInt_Check(_in)) return types::t_real(PyInt_AS_LONG(_in));
        if(PyFloat_Check(_in)) return types::t_real(PyFloat_AS_DOUBLE(_in));
        LADA_PYERROR(TypeError, "Expected quantity or number in input.");
        return types::t_real(0);
      }
      // creates global/local dictionary in order to run code.
      python::Object const globals = PyImport_ImportModule("quantities");
      if(not globals) return types::t_real(0);
      python::Object locals = PyDict_New();
      if(not locals) return types::t_real(0);
      python::Object units = PyString_FromString(_units.c_str());
      if(not units) return types::t_real(0);
      if(PyDict_SetItemString(locals.borrowed(), "number", _in) < 0)
        return types::t_real(0);
      if(PyDict_SetItemString(locals.borrowed(), "units", units.borrowed()) < 0)
        return types::t_real(0);
      python::Object const result(PyRun_String( "float(number.rescale(units))", Py_eval_input,
                                        PyModule_GetDict(globals.borrowed()), 
                                        locals.borrowed() ));
      if(not result) return types::t_real(0);
      types::t_real const realresult = PyFloat_AsDouble(result.borrowed());
      if(PyErr_Occurred()) return types::t_real(0);
      return realresult;
    }

  }
}
