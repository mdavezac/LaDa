//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#include "atat.hpp"
#include "atat.impl.hpp"

using namespace boost::python;

namespace LaDa
{
  namespace Python
  {
    void expose_atat() 
    {
      namespace bp = boost::python;
      expose_atatvector< atat::rVector3d >( "rVector3d", "a 3d-vector of real values." ); 
      expose_atatvector< atat::iVector3d >( "iVector3d", "a 3d-vector of integer values." ); 
      expose_atatmatrix< atat::rMatrix3d >( "rMatrix3d", "a 3x3 matrix of real values.\n" 
                                            " Note that the coefficients are accessed using a tuple"
                                            " as in \" a[(0,0)] \" where \"a\" is an rMatrix3d."); 
      bp::def( "inv_rMatrix3d",
               &LaDa::atat::details::inv_rMatrix3d< atat::rMatrix3d >,
               bp::arg("matrix"),
               "Inverts an rMatrix3d." );
      bp::def( "trans_rMatrix3d",
               &LaDa::atat::details::trans_rMatrix3d< atat::rMatrix3d >,
               bp::arg("matrix"),
               "Transpose of an rMatrix3d." );
    }
  }
} // namespace LaDa
