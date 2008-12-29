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
      expose_atatvector< atat::rVector3d >( "rVector3d" ); 
      expose_atatvector< atat::iVector3d >( "iVector3d" ); 
      expose_atatmatrix< atat::rMatrix3d >( "rMatrix3d" ); 
    }
  }
} // namespace LaDa
