#ifndef LADA_CRYSTAL_SCALE_H
#define LADA_CRYSTAL_SCALE_H

#include <Python.h>

#include <iostream>

#include <opt/types.h>

#include "object.h"
#include "exceptions.h"

namespace LaDa
{
  namespace python
  {
    //! \brief Wraps units from python package quantities.
    //! \throws Tries to only throw error::pyerror and derived.
    class Quantity : public Object
    {
        //! Private, can't happen cause it would have to throw. 
        Quantity(Object const &) {} 
      public:
        //! Acquires a python reference.
        Quantity(PyObject* _in = NULL) : Object(_in) {}
        //! \brief Sets quantity from given scale and units.
        //! \details Does not throw, but reference may be invalid after call.
        Quantity(types::t_real _in, std::string const &_units) : Object()
        {
          try { set(_in, _units); }
          catch(...) { object_ = NULL; }
        }
        //! Copy construct.
        Quantity(Quantity const &_in) : Object(_in) {}

        //! returns input in given units.
        types::t_real get(std::string const &_units = "")
        {
          Object const result
          (
            _units.empty() ? 
              run_code_("float(self.magnitude)"):
              run_code_("float(self.rescale('" + _units + "').magnitude)")
          );
          if(not result) BOOST_THROW_EXCEPTION(error::pyerror());
          return PyFloat_AS_DOUBLE(result.borrowed()); 
        }
        //! \brief Resets input to input scale.
        //! \details If input is a scale object, just change internal reference.
        //!          If it is an integer or float, then rescale while keeping to
        //!          the same units as currently held.
        void reset(Object const &_in)
        {
          if(not _in) this->Object::reset(NULL); 
          if(not Quantity::isinstance(_in)) this->Object::reset( run_code_("input*self.units", _in) );
          else this->Object::reset(_in); 
        }
  
        //! \brief True if input is an instance of
        //!        quantities.untiquantity.UnitQuantity or subtype.
        //! \details Never throws after python exception. Clears python error and
        //!          returns false.
        static bool isinstance(Object const &_in)
          { return Quantity::isinstance(_in.borrowed()); }
        //! \brief True if input is an instance of
        //!        quantities.untiquantity.UnitQuantity or subtype.
        //! \details Never throws after python exception. Clears python error and
        //!          returns false.
        static bool isinstance(PyObject *_in)
        {
          if(_in == NULL) return false;
          try
          { 
            return PyObject_IsInstance(_in, Quantity::classes_().borrowed()) == 1; 
          }
          catch(error::pyerror &_e) { PyErr_Clear(); return false; }
        }
  
        //! \brief New reference to quantities.unitquantity.Unitquantity class.
        //! \throws error::ImportError when the class cannot be imported.
        //!         Does not clear the pyerror exception. 
        static Object class_()
        {
          Object quant_type;
          Object quant( PyImport_ImportModule("quantities.unitquantity") );
          if(not quant) goto error; 
          quant_type.reset( PyObject_GetAttrString(quant.borrowed(), "UnitQuantity") );
          if(not quant_type) goto error; 
          return quant_type;
          
          error:
            PyErr_Clear();
            LADA_PYERROR(ImportError, "Could not import class Uniquantity from quantities package.");
            BOOST_THROW_EXCEPTION(error::ImportError());
        }
        //! \returns tuple containing UnitQuantity and Quantity classes.
        //! \throws error::ImportError when the class cannot be imported.
        //!         Does not clear the pyerror exception. 
        static Object classes_()
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
            LADA_PYERROR(ImportError, "Could not import class Uniquantity from quantities package.");
            BOOST_THROW_EXCEPTION(error::ImportError());
        }
        //! \brief New reference to quantity module.
        //! \throws error::ImportError when the package cannot be imported.
        //!         Does not clear the pyerror exception. 
        static Object module()
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

        //! Imports object from quantities package.
        static Object import(std::string const &_in)
        {
          Object const result = PyObject_GetAttrString(Quantity::module().borrowed(), _in.c_str()); 
          if(not result) BOOST_THROW_EXCEPTION(error::ImportError());
          return result;
        }
        //! Imports object from quantities package.
        static Object import(Object const &_in)
        {
          Object const result = PyObject_GetAttr(Quantity::module().borrowed(), _in.borrowed()); 
          if(not result) BOOST_THROW_EXCEPTION(error::ImportError());
          return result;
        }

      private:
        //! \brief Runs python code on current object
        //! \details The internal reference  is called self, whereas the input
        //!          object is called input if non null.
        Object run_code_(std::string const &_code, Object const &_in)
        {
          if(not is_valid())
          {
            LADA_PYERROR(RuntimeError, "Internal reference is not set.");
            BOOST_THROW_EXCEPTION(error::RuntimeError());
          }
          // creates global/local dictionary in order to run code.
          Object const globals = Quantity::module(); 
          Object locals = PyDict_New();
          if(not locals) BOOST_THROW_EXCEPTION(error::pyerror());
          if(_in and PyDict_SetItemString(locals.borrowed(), "input", _in.borrowed()) < 0)
            BOOST_THROW_EXCEPTION(error::pyerror()); 
          if(PyDict_SetItemString(locals.borrowed(), "self", object_) < 0)
            BOOST_THROW_EXCEPTION(error::pyerror());
          Object const result(PyRun_String( _code.c_str(), Py_eval_input,
                                            PyModule_GetDict(globals.borrowed()), 
                                            locals.borrowed() ));
          if(not result) BOOST_THROW_EXCEPTION(error::pyerror());
          return result;
        }
        //! \brief Runs python code on current object
        //! \details The internal reference  is called self.
        Object run_code_(std::string const &_code)
        {
          if(not is_valid())
          {
            LADA_PYERROR(RuntimeError, "Internal reference is not set.");
            BOOST_THROW_EXCEPTION(error::RuntimeError());
          }
          // creates global/local dictionary in order to run code.
          Object const globals = Quantity::module(); 
          Object locals = PyDict_New();
          if(not locals) BOOST_THROW_EXCEPTION(error::pyerror());
          if(PyDict_SetItemString(locals.borrowed(), "self", object_) < 0)
            BOOST_THROW_EXCEPTION(error::pyerror());
          Object const result(PyRun_String( _code.c_str(), Py_eval_input,
                                            PyModule_GetDict(globals.borrowed()), 
                                            locals.borrowed() ));
          if(not result) BOOST_THROW_EXCEPTION(error::pyerror());
          return result;
        }
        //! Sets internal object to represent the input double in given units.
        void set(double _in, std::string const &_units)
        {
          if(_units.size() == 0) BOOST_THROW_EXCEPTION(error::pyerror());
          // creates global/local dictionary in order to run code.
          Object globals = Quantity::module(); 
          Object locals = PyDict_New();
          if(not locals) BOOST_THROW_EXCEPTION(error::pyerror());
          Object infloat( PyFloat_FromDouble(_in) );
          if(not infloat) BOOST_THROW_EXCEPTION(error::pyerror());
          Object units(Quantity::import(_units));
          if(PyDict_SetItemString(locals.borrowed(), "scale", infloat.borrowed()) < 0)
            BOOST_THROW_EXCEPTION(error::pyerror());
          if(PyDict_SetItemString(locals.borrowed(), "units", units.borrowed()) < 0)
            BOOST_THROW_EXCEPTION(error::pyerror());
          Object const result( PyRun_String( "scale * units", Py_eval_input,
                                             PyModule_GetDict(globals.borrowed()), locals.borrowed()) );
          if(not result) BOOST_THROW_EXCEPTION(error::pyerror());
          this->Object::reset(result);
        }
    };
  }
}
#endif
