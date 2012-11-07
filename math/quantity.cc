#include "LaDaConfig.h"

#include "quantity.h"

namespace LaDa
{
  namespace math
  {
     types::t_real Quantity::get(std::string const &_units)
     {
       python::Object const result
       (
         _units.empty() ? 
           run_code_("float(self.magnitude)"):
           run_code_("float(self.rescale('" + _units + "').magnitude)")
       );
       if(not result) BOOST_THROW_EXCEPTION(error::internal());
       return PyFloat_AS_DOUBLE(result.borrowed()); 
     }

     void Quantity::reset(python::Object const &_in)
     {
       if(not _in) this->python::Object::reset(NULL); 
       if(not Quantity::isinstance(_in))
         this->python::Object::reset( run_code_("input*self.units", _in) );
       else this->python::Object::reset(_in); 
     }

     bool Quantity::isinstance(PyObject *_in)
     {
       if(_in == NULL) return false;
       try
       { 
         return PyObject_IsInstance(_in, Quantity::classes_().borrowed()) == 1; 
       }
       catch(error::root &_e) { PyErr_Clear(); return false; }
     }
     
     python::Object Quantity::class_()
     {
       python::Object quant_type;
       python::Object quant( PyImport_ImportModule("quantities.unitquantity") );
       if(not quant) goto error; 
       quant_type.reset( PyObject_GetAttrString(quant.borrowed(), "UnitQuantity") );
       if(not quant_type) goto error; 
       return quant_type;
       
       error:
         PyErr_Clear();
         LADA_PYTHROW(ImportError, "Could not import class Uniquantity from quantities package.");
     }

     python::Object Quantity::classes_()
     {
       python::Object quant_type, quant_type2;
       python::Object quant( PyImport_ImportModule("quantities.unitquantity") );
       if(not quant) goto error; 
       quant_type.reset( PyObject_GetAttrString(quant.borrowed(), "UnitQuantity") );
       if(not quant_type) goto error; 
       quant.reset(PyImport_ImportModule("quantities"));
       if(not quant) goto error;
       quant_type2.reset( PyObject_GetAttrString(quant.borrowed(), "Quantity") );
       if(not quant_type2) goto error;
       return PyTuple_Pack(2, quant_type2.borrowed(), quant_type.borrowed());

       error:
         PyErr_Clear();
         LADA_PYTHROW(ImportError, "Could not import class Uniquantity from quantities package.");
     }

     python::Object Quantity::module()
     {
       python::Object const quant( PyImport_ImportModule("quantities") );
       if(not quant)
       {
         PyErr_Clear();
         LADA_PYERROR(ImportError, "Could not import class Uniquantity from quantities package.");
         BOOST_THROW_EXCEPTION(error::ImportError());
       }
       return quant;
     }

     python::Object Quantity::import(std::string const &_in)
     {
       python::Object const result = PyObject_GetAttrString(Quantity::module().borrowed(), _in.c_str()); 
       if(not result) BOOST_THROW_EXCEPTION(error::ImportError());
       return result;
     }
     
     python::Object Quantity::import(python::Object const &_in)
     {
       python::Object const result = PyObject_GetAttr(Quantity::module().borrowed(), _in.borrowed()); 
       if(not result) BOOST_THROW_EXCEPTION(error::ImportError());
       return result;
     }

     python::Object Quantity::run_code_(std::string const &_code, python::Object const &_in)
     {
       if(not is_valid())
         LADA_PYTHROW(internal, "Internal reference is not set.");
       // creates global/local dictionary in order to run code.
       python::Object const globals = Quantity::module(); 
       python::Object locals = PyDict_New();
       if(not locals) BOOST_THROW_EXCEPTION(error::internal());
       if(_in and PyDict_SetItemString(locals.borrowed(), "input", _in.borrowed()) < 0)
         BOOST_THROW_EXCEPTION(error::internal()); 
       if(PyDict_SetItemString(locals.borrowed(), "self", object_) < 0)
         BOOST_THROW_EXCEPTION(error::internal());
       python::Object const result(PyRun_String( _code.c_str(), Py_eval_input,
                                         PyModule_GetDict(globals.borrowed()), 
                                         locals.borrowed() ));
       if(not result) BOOST_THROW_EXCEPTION(error::internal());
       return result;
     }
     //! \brief Runs python code on current object
     //! \details The internal reference  is called self.
     python::Object Quantity::run_code_(std::string const &_code)
     {
       if(not is_valid())
         LADA_PYTHROW(internal, "Internal reference is not set.");
       // creates global/local dictionary in order to run code.
       python::Object const globals = Quantity::module(); 
       python::Object locals = PyDict_New();
       if(not locals) BOOST_THROW_EXCEPTION(error::internal());
       if(PyDict_SetItemString(locals.borrowed(), "self", object_) < 0)
         BOOST_THROW_EXCEPTION(error::internal());
       python::Object const result(PyRun_String( _code.c_str(), Py_eval_input,
                                         PyModule_GetDict(globals.borrowed()), 
                                         locals.borrowed() ));
       if(not result) BOOST_THROW_EXCEPTION(error::internal());
       return result;
     }
     //! Sets internal object to represent the input double in given units.
     void Quantity::set(double _in, std::string const &_units)
     {
       if(_units.size() == 0) BOOST_THROW_EXCEPTION(error::internal());
       // creates global/local dictionary in order to run code.
       python::Object globals = Quantity::module(); 
       python::Object locals = PyDict_New();
       if(not locals) BOOST_THROW_EXCEPTION(error::internal());
       python::Object infloat( PyFloat_FromDouble(_in) );
       if(not infloat) BOOST_THROW_EXCEPTION(error::internal());
       python::Object units(Quantity::import(_units));
       if(PyDict_SetItemString(locals.borrowed(), "scale", infloat.borrowed()) < 0)
         BOOST_THROW_EXCEPTION(error::internal());
       if(PyDict_SetItemString(locals.borrowed(), "units", units.borrowed()) < 0)
         BOOST_THROW_EXCEPTION(error::internal());
       python::Object const result( PyRun_String( "scale * units", Py_eval_input,
                                          PyModule_GetDict(globals.borrowed()), locals.borrowed()) );
       if(not result) BOOST_THROW_EXCEPTION(error::internal());
       this->python::Object::reset(result);
     }

     types::t_real convert_toreal(PyObject *_in, std::string const &_units, types::t_real _default)
     {
       if(_in == NULL or _in == Py_None) return _default;
       return convert_toreal(_in, _units);
     }
     types::t_real convert_toreal(PyObject *_in, std::string const &_units)
     {
       if(PyObject_HasAttrString(_in, "rescale"))
       {
         char rescale_str[] = "rescale";
         char s_str[] = "s";
         python::Object rescaled 
             = PyObject_CallMethod( _in, rescale_str, s_str, _units.c_str());
         if(not rescaled)  return 0;
         python::Object magnitude = PyObject_GetAttrString(rescaled.borrowed(), "magnitude");
         if(not magnitude) return 0;
         char tolist_str[] = "tolist";
         python::Object real = PyObject_CallMethod(magnitude.borrowed(), tolist_str, NULL); 
         if(not real) return 0;
         _in = real.borrowed();
         if(PyInt_Check(_in) == 1)  return (types::t_real)PyInt_AS_LONG(_in);
         if(PyLong_Check(_in) == 1) return (types::t_real)PyLong_AsLong(_in);
         if(PyFloat_Check(_in) == 1) return PyFloat_AsDouble(_in);
       }
       else if(PyObject_HasAttrString(_in, "tolist"))
       {
         char tolist_str[] = "tolist";
         python::Object real = PyObject_CallMethod(_in, tolist_str, NULL); 
         if(not real) return 0;
         _in = real.borrowed();
       }
       if(PyInt_Check(_in) == 1)  return (types::t_real)PyInt_AS_LONG(_in);
       if(PyLong_Check(_in) == 1) return (types::t_real)PyLong_AsLong(_in);
       if(PyFloat_Check(_in) == 1) return PyFloat_AsDouble(_in);
       
       LADA_PYERROR(TypeError, "Could not convert value to real.");
       return 0;
     }
  }
}
