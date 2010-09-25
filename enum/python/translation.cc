#include "LaDaConfig.h"

#include <sstream>
#include <complex>

#include <boost/exception/diagnostic_information.hpp>

#include <boost/python/class.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/return_internal_reference.hpp>


#include <crystal/lattice.h>

#include "translation.hpp"
#include "../exceptions.h"
#include "../translation.h"



namespace LaDa
{
  
  namespace Python
  {
    namespace bp = boost::python;
    namespace e = enumeration;
    struct Translation : e::Translation
    {
      Translation   (math::iVector3d const &_smith, size_t _nsites)
                  : e::Translation(_smith, _nsites), first_(true) {}
      Translation   (const Translation &_c)
                  : e::Translation(_c), first_(_c.first_) {}

      Translation &iter()  { return *this; }
      Translation const& next()
      {
        if( first_ ) 
        {
          if( size() == 0 ) // no translations.
          {
            PyErr_SetString(PyExc_StopIteration, "End of translations.");
            bp::throw_error_already_set();
          }
          else first_ = false; 
        }
        else if( not e::Translation::operator++() )
        {
          first_ = true;
          PyErr_SetString(PyExc_StopIteration, "End of translations.");
          bp::throw_error_already_set();
        }
        return *this;
      }
      bool first_;
    };

    e::t_uint __call__(Translation const &_self, e::t_uint x, e::FlavorBase const &_fl) 
    {
      try { return _self(x, _fl); }
      catch( e::internal & x )
      {
        if( std::string const * mi=boost::get_error_info<e::error_string>(x) )
          PyErr_SetString(PyExc_RuntimeError, ("Error in translation: " + (*mi)).c_str());
        else
          PyErr_SetString(PyExc_RuntimeError, "Error in translation.");
      }
      catch( e::integer_too_large & x )
      {
        if( std::string const * mi=boost::get_error_info<e::error_string>(x) )
          PyErr_SetString( PyExc_ValueError, 
                           ("Error in translation: integer too large.\n"+(*mi)).c_str() );
        else
          PyErr_SetString(PyExc_ValueError, "Error in translation: integer too large.\n");
      }
      catch( e::argument_error & x )
      {
        if( std::string const * mi=boost::get_error_info<e::error_string>(x) )
          PyErr_SetString(PyExc_ValueError, ("Error in translation: " + (*mi)).c_str());
        else
          PyErr_SetString(PyExc_ValueError, "Error in translation.");
      }
      bp::throw_error_already_set();
      return -1;
    }

    void expose_translation()
    {
      bp::class_<Translation>
      (
        "Translation", 
        "Rotation + translation.", 
        bp::init<math::iVector3d const &, size_t>()
      ).def( bp::init<Translation const&>() )
       .def("__call__", &__call__)
       .def("__len__", &Translation::size)
       .def("__iter__", &Translation::iter, bp::return_internal_reference<1>())
       .def("next", &Translation::next, bp::return_internal_reference<1>());
    }

  }
} // namespace LaDa
