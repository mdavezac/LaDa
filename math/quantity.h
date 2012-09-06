#ifndef LADA_CRYSTAL_QUANTITY_H
#define LADA_CRYSTAL_QUANTITY_H

#include "LaDaConfig.h"
#include <Python.h>

#include <iostream>

#include <boost/python/errors.hpp>
#include <opt/types.h>

#include <python/object.h>
#include <python/exceptions.h>

namespace LaDa
{
  namespace math
  {
    //! \brief Wraps units from python package quantities.
    class Quantity : public python::Object
    {
        //! Private, can't happen cause it would have to throw. 
        Quantity(python::Object const &) {} 
      public:
        //! Acquires a python reference.
        Quantity(PyObject* _in = NULL) : python::Object(_in) {}
        //! \brief Sets quantity from given scale and units.
        //! \details Does not throw, but reference may be invalid after call.
        Quantity(types::t_real _in, std::string const &_units) : python::Object()
        {
          try { set(_in, _units); }
          catch(...) { object_ = NULL; }
        }
        //! Copy construct.
        Quantity(Quantity const &_in) : python::Object(_in) {}

        //! returns input in given units.
        types::t_real get(std::string const &_units = "");
        //! \brief Resets input to input scale.
        //! \details If input is a scale object, just change internal reference.
        //!          If it is an integer or float, then rescale while keeping to
        //!          the same units as currently held.
        void reset(python::Object const &_in);
  
        //! \brief True if input is an instance of
        //!        quantities.untiquantity.UnitQuantity or subtype.
        //! \details Never throws after python exception. Clears python error and
        //!          returns false.
        static bool isinstance(python::Object const &_in)
          { return Quantity::isinstance(_in.borrowed()); }
        //! \brief True if input is an instance of
        //!        quantities.untiquantity.UnitQuantity or subtype.
        //! \details Never throws after python exception. Clears python error and
        //!          returns false.
        static bool isinstance(PyObject *_in);
  
        //! \brief New reference to quantities.unitquantity.Unitquantity class.
        //! \throws error::ImportError when the class cannot be imported.
        //!         Does not clear the exception. 
        static python::Object class_();
        //! \returns tuple containing UnitQuantity and Quantity classes.
        //! \throws error::ImportError when the class cannot be imported.
        //!         Does not clear the exception. 
        static python::Object classes_();
        //! \brief New reference to quantity module.
        //! \throws error::ImportError when the package cannot be imported.
        //!         Does not clear the  exception. 
        static python::Object module();

        //! Imports object from quantities package.
        static python::Object import(std::string const &_in);
        //! Imports object from quantities package.
        static python::Object import(python::Object const &_in);

      private:
        //! \brief Runs python code on current object
        //! \details The internal reference  is called self, whereas the input
        //!          object is called input if non null.
        python::Object run_code_(std::string const &_code, python::Object const &_in);
        //! \brief Runs python code on current object
        //! \details The internal reference  is called self.
        python::Object run_code_(std::string const &_code);
        //! Sets internal object to represent the input double in given units.
        void set(double _in, std::string const &_units);
    };
    //! Converts an object to a real number. 
    //! \param[in] _in: Can be a float, a long, a int, a scalar numpy array, or
    //!                 a signed scalar array.
    //! \param[in] _units: The units to which a signed array should be converted.
    //! \param[in] _defaults: A default value if None or NULL is given.
    types::t_real convert_toreal(PyObject *_in, std::string const &_units, types::t_real _default);
  }
}
#endif
