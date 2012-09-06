#include "LaDaConfig.h"

#include "quantity.h"

namespace LaDa
{
  namespace python
  {
     types::t_real Quantity::get(std::string const &_units = "")
     {
       Object const result
       (
         _units.empty() ? 
           run_code_("float(self.magnitude)"):
           run_code_("float(self.rescale('" + _units + "').magnitude)")
       );
       if(not result) BOOST_THROW_EXCEPTION(error::internal());
       return PyFloat_AS_DOUBLE(result.borrowed()); 
     }

     void Quantity::reset(Object const &_in)
     {
       if(not _in) this->Object::reset(NULL); 
       if(not Quantity::isinstance(_in))
         this->Object::reset( run_code_("input*self.units", _in) );
       else this->Object::reset(_in); 
     }

     static bool Quantity::isinstance(PyObject *_in)
     {
       if(_in == NULL) return false;
       try
       { 
         return PyObject_IsInstance(_in, Quantity::classes_().borrowed()) == 1; 
       }
       catch(error::root &_e) { PyErr_Clear(); return false; }
     }
     
     static Object Quantity::class_()
     {
       Object quant_type;
       Object quant( PyImport_ImportModule("quantities.unitquantity") );
       if(not quant) goto error; 
       quant_type.reset( PyObject_GetAttrString(quant.borrowed(), "UnitQuantity") );
       if(not quant_type) goto error; 
       return quant_type;
       
       error:
         PyErr_Clear();
         LADA_PYTHROW(ImportError, "Could not import class Uniquantity from quantities package.");
     }

     static Object Quantity::classes_()
     {
       Object quant_type, quant_type2;
       Object quant( PyImport_ImportModule("quantities.unitquantity") );
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

     static Object Quantity::module()
     {
       Object const quant( PyImport_ImportModule("quantities") );
       if(not quant)
       {
         PyErr_Clear();
         LADA_PYERROR(ImportError, "Could not import class Uniquantity from quantities package.");
         BOOST_THROW_EXCEPTION(error::ImportError());
       }
       return quant;
     }

     static Object Quantity::import(std::string const &_in)
     {
       Object const result = PyObject_GetAttrString(Quantity::module().borrowed(), _in.c_str()); 
       if(not result) BOOST_THROW_EXCEPTION(error::ImportError());
       return result;
     }
     
     static Object Quantity::import(Object const &_in)
     {
       Object const result = PyObject_GetAttr(Quantity::module().borrowed(), _in.borrowed()); 
       if(not result) BOOST_THROW_EXCEPTION(error::ImportError());
       return result;
     }

     Object Quantity::run_code_(std::string const &_code, Object const &_in)
     {
       if(not is_valid())
         LADA_PYTHROW(internal, "Internal reference is not set.");
       // creates global/local dictionary in order to run code.
       Object const globals = Quantity::module(); 
       Object locals = PyDict_New();
       if(not locals) BOOST_THROW_EXCEPTION(error::internal());
       if(_in and PyDict_SetItemString(locals.borrowed(), "input", _in.borrowed()) < 0)
         BOOST_THROW_EXCEPTION(error::internal()); 
       if(PyDict_SetItemString(locals.borrowed(), "self", object_) < 0)
         BOOST_THROW_EXCEPTION(error::internal());
       Object const result(PyRun_String( _code.c_str(), Py_eval_input,
                                         PyModule_GetDict(globals.borrowed()), 
                                         locals.borrowed() ));
       if(not result) BOOST_THROW_EXCEPTION(error::internal());
       return result;
     }
     //! \brief Runs python code on current object
     //! \details The internal reference  is called self.
     Object Quantity::run_code_(std::string const &_code)
     {
       if(not is_valid())
         LADA_PYTHROW(internal, "Internal reference is not set.");
       // creates global/local dictionary in order to run code.
       Object const globals = Quantity::module(); 
       Object locals = PyDict_New();
       if(not locals) BOOST_THROW_EXCEPTION(error::internal());
       if(PyDict_SetItemString(locals.borrowed(), "self", object_) < 0)
         BOOST_THROW_EXCEPTION(error::internal());
       Object const result(PyRun_String( _code.c_str(), Py_eval_input,
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
       Object globals = Quantity::module(); 
       Object locals = PyDict_New();
       if(not locals) BOOST_THROW_EXCEPTION(error::internal());
       Object infloat( PyFloat_FromDouble(_in) );
       if(not infloat) BOOST_THROW_EXCEPTION(error::internal());
       Object units(Quantity::import(_units));
       if(PyDict_SetItemString(locals.borrowed(), "scale", infloat.borrowed()) < 0)
         BOOST_THROW_EXCEPTION(error::internal());
       if(PyDict_SetItemString(locals.borrowed(), "units", units.borrowed()) < 0)
         BOOST_THROW_EXCEPTION(error::internal());
       Object const result( PyRun_String( "scale * units", Py_eval_input,
                                          PyModule_GetDict(globals.borrowed()), locals.borrowed()) );
       if(not result) BOOST_THROW_EXCEPTION(error::internal());
       this->Object::reset(result);
     }

     types::t_real convert_toreal(PyObject *_in, std::string const &_units, types::t_real _default)
     {
       if(PyObject_HasAttrString(_in, "rescale"))
       {
         python::Object rescaled 
             = PyObject_CallMethod(_in, "rescale", "s", _units.c_str());
         if(not rescaled)  return NULL;
         python::Object magnitude = PyObject_GetAttrString(rescaled.borrowed(), "magnitude");
         if(not magnitude) return NULL;
         python::Object real = PyObject_CallMethodObjArgs(magnitude.borrowed(), "tolist", NULL);
         if(not real) return NULL;
         _in = real.borrowed();
       }
       else if(_in == NULL or _in == Py_None) return _default;
       if(PyInt_Check(_in))  return (types::t_real)PyInt_AS_LONG(_in);
       if(PyLong_Check(_in)) return (types::t_real)PyLong_AS_LONG(_in);
       
       LADA_PYERROR(TypeError, "cutoff should be a floating point or an energy.");
       return -1;
     }
  }
}
