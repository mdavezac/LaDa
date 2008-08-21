
//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <sstream>

#include <boost/python.hpp>

#include <opt/types.h>
#include <opt/errors.h>

#include "errortuple.hpp"


namespace PythonLaDa
{
  std::string print_error( const opt::ErrorTuple &_error )
  {
    std::ostringstream sstr;
    sstr << _error;
    return sstr.str();
  }
  opt::ErrorTuple* add( const opt::ErrorTuple &_a, const opt::ErrorTuple &_b )
  {
    opt::ErrorTuple *result = new opt::ErrorTuple( _a );
    *result += _b;
    return result;
  }
  void expose_errors()
  {
    using namespace boost::python;
    class_< opt::ErrorTuple >( "ErrorTuple" )
      .def( init< opt::ErrorTuple >() )
      .def( init< types::t_real, types::t_real >() )
      .def( "mean", &opt::ErrorTuple::mean )
      .def( "variance", &opt::ErrorTuple::variance )
      .def( "max", &opt::ErrorTuple::max )
      .def( self += other< opt::ErrorTuple >() ) 
      .def( "__add__", &add,
            return_value_policy< manage_new_object >() )
      .def( "__str__",  &print_error );
  }
}
