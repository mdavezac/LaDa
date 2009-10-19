//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <iostream>
#include <sstream>
#include <complex>

#include <boost/python/def.hpp>
#include <pyublas/numpy.hpp>


#include "find_pis.hpp"
#include "../find_pis.h"


namespace LaDa
{
  namespace Python
  {
    pyublas::numpy_vector<types::t_real> find_pis( CE::t_ClusterClasses const &_cls,
                                                   Crystal::Structure const &_str,
                                                   size_t _site )
    {
      pyublas::numpy_vector<types::t_real> vec; 
      CE::find_pis(_cls, _str, vec, _site);
      return vec;
    }

    void expose_find_pis()
    {
      namespace bp = boost::python;

      bp::def
      (
        "find_pis",
        &find_pis,
        ( bp::arg("classes"), bp::arg("structure"), bp::arg("site")=0 ),
        "Returns pis of a given structure for a given array of classes of equivalent figures,\n"
        "with site the site index of the origin of the figure."
      );
    }

  }
} // namespace LaDa
