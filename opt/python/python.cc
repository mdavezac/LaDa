//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>

#include <revision.h>

#include "convexhull.impl.hpp"
#include "errortuple.hpp"

BOOST_PYTHON_MODULE(opt)
{
  LaDa::Python::expose_errors();
  LaDa::Python::exposeConvexHull<boost::python::object>( "ConvexHull" );
}
