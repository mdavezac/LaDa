#include "LaDaConfig.h"

#include <boost/python.hpp>
#include <boost/python/module.hpp>

#include "clj.hpp"
#include "functional.hpp"

BOOST_PYTHON_MODULE(models)
{
  LaDa::Python::expose_clj();
  LaDa::Python::expose_functional();
}
