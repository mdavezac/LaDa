#ifndef LADA_CRYSTAL_SCALE_H
#define LADA_CRYSTAL_SCALE_H

#include <Python.h>

#include <iostream>

#include <boost/python/errors.hpp>
#include <opt/types.h>

#include "object.h"
#include "exceptions.h"

namespace LaDa
{
  namespace python
  {
    //! \brief Wraps units from python package quantities.
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
        types::t_real get(std::string const &_units = "");
        //! \brief Resets input to input scale.
        //! \details If input is a scale object, just change internal reference.
        //!          If it is an integer or float, then rescale while keeping to
        //!          the same units as currently held.
        void reset(Object const &_in);
  
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
        static bool isinstance(PyObject *_in);
  
        //! \brief New reference to quantities.unitquantity.Unitquantity class.
        //! \throws error::ImportError when the class cannot be imported.
        //!         Does not clear the exception. 
        static Object class_();
        //! \returns tuple containing UnitQuantity and Quantity classes.
        //! \throws error::ImportError when the class cannot be imported.
        //!         Does not clear the exception. 
        static Object classes_();
        //! \brief New reference to quantity module.
        //! \throws error::ImportError when the package cannot be imported.
        //!         Does not clear the  exception. 
        static Object module();

        //! Imports object from quantities package.
        static Object import(std::string const &_in);
        //! Imports object from quantities package.
        static Object import(Object const &_in);

      private:
        //! \brief Runs python code on current object
        //! \details The internal reference  is called self, whereas the input
        //!          object is called input if non null.
        Object run_code_(std::string const &_code, Object const &_in);

        //! \brief Runs python code on current object
        //! \details The internal reference  is called self.
        Object run_code_(std::string const &_code);
        
        //! Sets internal object to represent the input double in given units.
        void set(double _in, std::string const &_units);
    };

    //! Converts python object to real.
    //! \param[in] _in: Python object which can be a scalar array, a float, or
    //!                 an integer. It can also be signed by a quantity.
    //! \param[in] _units: Units to which the value should be converted if signed by a unit.                
    //! \param[in] _default: Default value if NULL or None.
    types::t_real convert_toreal(PyObject *_in, std::string const &_units, types::t_real _default);
  }
}
#endif
