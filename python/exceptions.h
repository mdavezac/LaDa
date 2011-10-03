#ifndef LADA_PYTHON_EXCEPTIONS_H
#define LADA_PYTHON_EXCEPTIONS_H
#include "LaDaConfig.h"

#include <boost/python/object.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/str.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/exception_translator.hpp>
#include <boost/exception/diagnostic_information.hpp>

#include <root_exceptions.h>

namespace LaDa
{
  namespace error
  {
    //! Attribute error thrown explicitely by lada.
    struct AttributeError: virtual input {};
    //! Key error thrown explicitely by lada.
    struct KeyError: virtual input, virtual out_of_range {};
    //! Value error thrown explicitely by lada.
    struct ValueError: virtual input {};
    //! Index error thrown explicitely by lada.
    struct IndexError: virtual input {};
    //! Argument error thrown explicitely by lada.
    struct TypeError: virtual input {};
    //! Not implemented error thrown explicitely by lada.
    struct NotImplementedError: virtual internal {};
    //! Internal error thrown explicitely by lada.
    struct InternalError: virtual internal {};
  }

  namespace python
  {
    // Class to declare and register c++ to python exceptions.
    template<class T> 
      class PyException
      {
        public:
          PyException() {}
          boost::python::object const &initialize( std::string const &_name, 
                                                   std::string const &_doc, 
                 boost::python::tuple const &_bases
                   = boost::python::make_tuple(boost::python::object(
                       boost::python::handle<>(
                         PyExc_StandardError))) ) 
          {
            static bool is_first = true;
            if(is_first)
            {
              namespace bp = boost::python;
              bp::dict d;
              d["__doc__"] = bp::str(_doc);
              d["name"] = bp::str(_name);
              name_ = _name;
              doc_ = _doc;
              exception_ = bp::object(bp::handle<>(bp::borrowed(
                PyErr_NewExceptionWithDoc(&name_[0], &doc_[0], _bases.ptr(), d.ptr()) ) ) );
              is_first = false;
            }
            return exception_;
          }

          void operator()(T const &_e) const
          {
            std::string message = boost::diagnostic_information(_e);
            if(    name_[0] == 'a' or name_[0] == 'e' or name_[0] == 'i'
                or name_[0] == 'o' or name_[0] == 'u' or name_[0] == 'y' )
              message += "Encountered an " + name_ + " error.";
            else 
              message += "Encountered a " + name_ + " error.";
            PyErr_SetString(exception_.ptr(), message.c_str());
          }

          static void throw_error(std::string const &_message)
          {
            PyErr_SetString(exception_.ptr(), _message.c_str());
            boost::python::throw_error_already_set();
          };

          static boost::python::object const & exception() { return exception_; }

        private:
          std::string name_;
          std::string doc_;
          //! static exception object. Act as global.
          static boost::python::object exception_;
      };
    template<class T> boost::python::object PyException<T>::exception_ = boost::python::object();
#   ifdef LADA_REGISTER_PYEXCEPT
#     error LADA_REGISTER_PYEXCEPT already defined.
#   endif
#   ifdef LADA_REGISTER_PYEXCEPT_WITH_BASE
#     error LADA_REGISTER_PYEXCEPT_WITH_BASE already defined.
#   endif
#   define LADA_REGISTER_PYEXCEPT(T, n, doc, scope)\
    { \
      std::string name = n; \
      ::LaDa::python::PyException<T> e; \
      e.initialize(name, doc); \
      boost::python::register_exception_translator<T>(e);\
      scope.attr(name.substr(name.rfind('.')+1).c_str()) = ::LaDa::python::PyException<T>::exception();\
    }
#   define LADA_REGISTER_PYEXCEPT_WITH_BASE(T, n, doc, scope, base)\
    { \
      std::string name = n; \
      ::LaDa::python::PyException<T> e; \
      e.initialize(name, doc, base); \
      boost::python::register_exception_translator<T>(e);\
      scope.attr(name.substr(name.rfind('.')+1).c_str()) = ::LaDa::python::PyException<T>::exception();\
    }
  }
}
# endif 
