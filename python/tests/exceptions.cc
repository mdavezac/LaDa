#include "LaDaConfig.h"

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

#include <root_exceptions.h>
#include "../exceptions.h"

namespace bp = boost::python;

void nothrow () {};
void dothrow_nomessage() { BOOST_THROW_EXCEPTION(::LaDa::error::LADA_TYPE()); }
void dothrow_message()
{
  std::string message = "This is a message.";
  BOOST_THROW_EXCEPTION(::LaDa::error::LADA_TYPE() << LaDa::error::string(message)); 
}


BOOST_PYTHON_MODULE(LADA_MODULE)
{
  bp::def("nothrow", &nothrow);
  bp::def("dothrow_nomessage", &dothrow_nomessage);
  bp::def("dothrow_message", &dothrow_message);
}
