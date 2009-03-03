//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>

#include <revision.h>

#include <python/std_vector.hpp>
#include <opt/types.h>

#include "convexhull.impl.hpp"
#include "errortuple.hpp"

BOOST_PYTHON_MODULE(opt)
{
  LaDa::Python::expose_errors();
  LaDa::Python::exposeConvexHull<boost::python::object>( "ConvexHull" );
  LaDa::Python::expose_vector<LaDa::types::t_real>( "cReals", "A stl container of real values." );
}
