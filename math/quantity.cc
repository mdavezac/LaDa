#include "LaDaConfig.h"

#include <errors/exceptions.h>
#include "quantity.h"

#include <iostream>


namespace LaDa
{
  namespace math
  {
    static PyObject* UnitQuantityClass()
    {
      PyObject* quant( PyImport_ImportModule("quantities") );
      if(not quant) return NULL;
      PyObject* result = PyObject_GetAttrString(quant, "UnitQuantity");
      Py_DECREF(quant);
      return result;
    }
    static PyObject* QuantityClass()
    {
      PyObject* quant( PyImport_ImportModule("quantities") );
      if(not quant) return NULL;
      PyObject* result = PyObject_GetAttrString(quant, "Quantity");
      Py_DECREF(quant);
      return result;
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
      PyObject* globals = PyImport_ImportModule("quantities");
      if(not globals) return NULL;
      PyObject* locals = PyDict_New();
      if(not locals) {Py_DECREF(globals); return NULL; }
      PyObject* number = PyFloat_FromDouble(_double);
      if(not number) {Py_DECREF(globals); Py_DECREF(number); return NULL;}
      PyObject* units = PyString_FromString(_units.c_str());
      if(not units) return NULL;
      if(PyDict_SetItemString(locals, "number", number) < 0)
        {Py_DECREF(globals); Py_DECREF(locals); Py_DECREF(number); Py_DECREF(units); return NULL;}
      Py_DECREF(number);
      if(PyDict_SetItemString(locals, "units", units) < 0)
        {Py_DECREF(globals); Py_DECREF(locals); Py_DECREF(units); return NULL;}
      Py_DECREF(units);
      PyObject* result(PyRun_String( "quantity.Quantity(number, units)", Py_eval_input,
                                     PyModule_GetDict(globals), 
                                     locals ));
      Py_DECREF(locals);
      Py_DECREF(globals);
      return result;
    }

    PyObject *PyQuantity_FromCWithTemplate(types::t_real const &_double, PyObject *_unittemplate)
    {
      // creates global/local dictionary in order to run code.
      PyObject* globals = PyImport_ImportModule("quantities");
      if(not globals) return NULL;
      PyObject* locals = PyDict_New();
      if(not locals) {Py_DECREF(globals); return NULL; }
      PyObject* number = PyFloat_FromDouble(_double);
      if(not number) {Py_DECREF(globals); Py_DECREF(number); return NULL;}
      if(PyDict_SetItemString(locals, "number", number) < 0)
        {Py_DECREF(globals); Py_DECREF(locals); Py_DECREF(number); return NULL;}
      Py_DECREF(number);
      if(PyDict_SetItemString(locals, "units", _unittemplate) < 0)
        {Py_DECREF(globals); Py_DECREF(locals); return NULL;}
      PyObject* result(PyRun_String( "quantity.Quantity(number, units.units)", Py_eval_input,
                                     PyModule_GetDict(globals), 
                                     locals ));
      Py_DECREF(locals);
      Py_DECREF(globals);
      return result;
    }

    static bool PyQuantity_Convertible(PyObject *_a, PyObject *_b)
    {
      if(not PyQuantity_Check(_a)) return false;
      if(not PyQuantity_Check(_b)) return false;
      char rescale_str[] = "rescale";
      char s_str[] = "O";
      PyObject* units = PyObject_GetAttrString(_b, "units");
      if(not units) return false;
      PyObject *result = PyObject_CallMethod(_a, rescale_str, s_str, units);
      Py_DECREF(units);
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
        return _number;
      }
      
      // creates global/local dictionary in order to run code.
      PyObject* locals = PyDict_New();
      if(not locals) return NULL;
      if(PyDict_SetItemString(locals, "number", _number) < 0)
        { Py_DECREF(locals); return NULL; }
      if(PyDict_SetItemString(locals, "unitdims", _units) < 0)
        { Py_DECREF(locals); return NULL; }
      PyObject* globals = PyImport_ImportModule("quantities");
      if(not globals) {Py_DECREF(locals); return NULL;}
      PyObject* result = PyRun_String( "Quantity(number, unitdims.units)",
                                       Py_eval_input,
                                       PyModule_GetDict(globals), 
                                       locals );
      Py_DECREF(locals);
      Py_DECREF(globals);
      return result;
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
      PyObject* globals = PyImport_ImportModule("quantities");
      if(not globals) return types::t_real(0);
      PyObject* locals = PyDict_New();
      if(not locals) {Py_DECREF(globals); return types::t_real(0);}
      PyObject* units = PyString_FromString(_units.c_str());
      if(not units) {Py_DECREF(locals); Py_DECREF(globals); return types::t_real(0); }
      if(PyDict_SetItemString(locals, "units", units) < 0)
        {Py_DECREF(locals); Py_DECREF(globals); Py_DECREF(units); return types::t_real(0); }
      Py_DECREF(units);
      if(PyDict_SetItemString(locals, "number", _in) < 0)
        {Py_DECREF(locals); Py_DECREF(globals); return types::t_real(0); }
      PyObject* result = PyRun_String( "float(number.rescale(units))", Py_eval_input,
                                       PyModule_GetDict(globals), 
                                       locals );
      Py_DECREF(globals);
      Py_DECREF(locals);
      if(not result) return types::t_real(0);
      types::t_real const realresult = PyFloat_AsDouble(result);
      Py_DECREF(result);
      if(PyErr_Occurred()) return types::t_real(0);
      return realresult;
    }
    types::t_real PyQuantity_GetPy(PyObject *_number, PyObject *_units)
    {
      if(not PyQuantity_Check(_number))
      {
        LADA_PYERROR(TypeError, "PyQuantity_GetPy: First argument should be a quantity.");
        return types::t_real(0);
      }
      if(not PyQuantity_Check(_units))
      {
        LADA_PYERROR(TypeError, "PyQuantityGetPy: Second argument should be a quantity.");
        return types::t_real(0);
      }
      // creates global/local dictionary in order to run code.
      PyObject* locals = PyDict_New();
      if(not locals) return types::t_real(0);
      if(PyDict_SetItemString(locals, "number", _number) < 0)
        {Py_DECREF(locals); return types::t_real(0);}
      if(PyDict_SetItemString(locals, "units", _units) < 0)
        {Py_DECREF(locals); return types::t_real(0);}
      PyObject* globals = PyImport_ImportModule("quantities");
      if(not globals) {Py_DECREF(locals); return types::t_real(0);}
      PyObject* result = PyRun_String( "float(number.rescale(units.units))",
                                       Py_eval_input,
                                       PyModule_GetDict(globals), 
                                       locals );
      Py_DECREF(globals);
      Py_DECREF(locals);
      if(not result) return types::t_real(0);
      types::t_real const realresult = PyFloat_AsDouble(result);
      Py_DECREF(result);
      if(PyErr_Occurred()) return types::t_real(0);
      return realresult;
    }

    types::t_real PyQuantity_AsReal(PyObject *_in)
    {
      PyObject* locals = PyDict_New();
      if(not locals) return types::t_real(0);
      if(PyDict_SetItemString(locals, "number", _in) < 0)
        { Py_DECREF(locals); return types::t_real(0);}
      PyObject *globals = PyEval_GetBuiltins();
      if(not globals) {Py_DECREF(locals); return types::t_real(0);}
      PyObject* result = PyRun_String( "float(number)", Py_eval_input,
                                       globals, locals );
      Py_DECREF(locals);
      Py_DECREF(globals);
      if(not result) return types::t_real(0);
      types::t_real const realresult = PyFloat_AsDouble(result);
      Py_DECREF(result);
      if(PyErr_Occurred()) return types::t_real(0);
      return realresult;
    }
  }
}
