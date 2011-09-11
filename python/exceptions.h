#ifndef LADA_PYTHON_EXCEPTIONS_H
#define LADA_PYTHON_EXCEPTIONS_H
#include "LaDaConfig.h"

#include <boost/python/object.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/str.hpp>
#include <boost/python/exception_translator.hpp>

namespace LaDa
{
  namespace python
  {
    // Class to declare and register c++ to python exceptions.
    template<class T> 
      class PyException
      {
        public:
          PyException   (std::string const &_name, std::string _doc)
                      : name_(_name), doc_(_doc) {}; 
          boost::python::object const &get() const
          {
            static std::string name = name_;
            static std::string doc = doc_;
            static boost::python::object exception_(boost::python::handle<>(
              PyErr_NewExceptionWithDoc(&name[0], &doc[0], NULL, NULL) ) );
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
            PyErr_SetString(get().ptr(), message.c_str());
          }

        private:
          std::string name_;
          std::string doc_;
      };
#   ifdef LADA_REGISTER_PYEXCEPT
#     error LADA_REGISTER_PYEXCEPT already defined.
#   endif
#   define LADA_REGISTER_PYEXCEPT(T, n, doc, scope)\
    { \
      std::string name = n; \
      ::LaDa::python::PyException<T> e(name, doc); \
      boost::python::register_exception_translator<T>(e);\
      if(name.rfind('.') != std::string::npos)\
        scope.attr(name.substr(name.rfind('.')+1).c_str()) = e.get();\
      else scope(name) = e.get();\
    }
  }
}
# endif 
