#include "PyladaConfig.h"

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

#include <root_exceptions.h>
#include "../exceptions.h"

namespace bp = boost::python;

void nothrow () {};
void dothrow_nomessage() {
  BOOST_THROW_EXCEPTION(::Pylada::error::PYLADA_TYPE()); 
}
void dothrow_message()
{
  std::string message = "This is a message.";
  BOOST_THROW_EXCEPTION(::Pylada::error::PYLADA_TYPE() << Pylada::error::string(message)); 
}
void dopythrow_message()
{
  PyErr_SetString(::Pylada::error::get_error(PYLADA_TYPENAME).ptr(), "This is another message.");
  std::string message = "This is a message.";
  BOOST_THROW_EXCEPTION(::Pylada::error::PYLADA_TYPE() << Pylada::error::string(message)); 
}

//void BOOST_THROW_EXCEPTION( error::root exc) {
//  throw exc.what();
//}

BOOST_PYTHON_MODULE(PYLADA_MODULE)
{
  Pylada::error::bp_register();
  bp::def("nothrow", &nothrow);
  bp::def("dothrow_nomessage", &dothrow_nomessage);
  bp::def("dothrow_message", &dothrow_message);
  bp::def("dopythrow_message", &dopythrow_message);
}
