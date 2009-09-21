//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#include <boost/python/module_init.hpp>

#include "atat.impl.hpp"

#include "../is_int.h"
using namespace boost::python;

namespace LaDa
{
  namespace Python
  {
    template< class T_MATRIX > types::t_real det( const T_MATRIX& _a )
      { return atat::det( _a ); }

    bool is_vinteger( atat::rVector3d const &_a ) 
      {return atat::is_integer(_a); }
    bool is_minteger( atat::rMatrix3d const &_a ) 
      {return atat::is_integer(_a); }

    void expose_atat() 
    {
      namespace bp = boost::python;
      expose_atatvector< atat::rVector3d >( "rVector3d", "a 3d-vector of real values." ); 
      expose_atatvector< atat::iVector3d >( "iVector3d", "a 3d-vector of integer values." ); 
      expose_atatmatrix< atat::rMatrix3d >( "rMatrix3d", "a 3x3 matrix of real values.\n" 
                                            " Note that the coefficients are accessed using a tuple"
                                            " as in \" a[(0,0)] \" where \"a\" is an rMatrix3d."); 
      bp::def( "inverse",
               &LaDa::atat::details::inv_rMatrix3d< atat::rMatrix3d >,
               bp::arg("matrix"),
               "Inverts an rMatrix3d." );
      bp::def( "transpose",
               &LaDa::atat::details::trans_rMatrix3d< atat::rMatrix3d >,
               bp::arg("matrix"),
               "Transpose of an rMatrix3d." );
      bp::def( "det", &det< atat::rMatrix3d >, bp::arg("matrix"), "Returns determinant." );

      boost::python::def("is_integer", &is_vinteger);
      boost::python::def("is_integer", &is_minteger);
    }

  }
} // namespace LaDa
BOOST_PYTHON_MODULE(atat)
{
  // PythonLaDa::expose_svn();
  LaDa::Python::expose_atat();
}
