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
#include <boost/python/list.hpp>
#include <pyublas/numpy.hpp>


#include "find_pis.hpp"
#include "../find_pis.h"


namespace LaDa
{
  namespace Python
  {
    namespace bp = boost::python;
    pyublas::numpy_vector<types::t_real> find_pis( CE::t_ClusterClasses const &_cls,
                                                   Crystal::Structure const &_str,
                                                   size_t _site )
    {
      pyublas::numpy_vector<types::t_real> vec; 
      CE::find_pis(_cls, _str, vec, _site);
      return vec;
    }
    pyublas::numpy_vector<types::t_real> find_pis3( CE::t_ClusterClasses const &_cls,
                                                   Crystal::Structure const &_str)
    {
      pyublas::numpy_vector<types::t_real> vec; 
      CE::find_pis2(_cls, _str, vec);
      return vec;
    }

    pyublas::numpy_vector<types::t_real> find_pis2( bp::list const &_cls,
                                                    Crystal::Structure const &_str,
                                                    size_t _site )
    {
      CE::t_ClusterClasses cls(bp::len(_cls));
      for(size_t i(0); i < bp::len(_cls); ++i)
        cls[i] = bp::extract<CE::t_Clusters>(_cls[i]);
      return find_pis(cls, _str, _site);
    }

    void expose_find_pis()
    {

      bp::def
      (
        "find_pis",
        &find_pis,
        ( bp::arg("classes"), bp::arg("structure"), bp::arg("site")=0 ),
        "Returns pis of a given structure for a given array of classes of equivalent figures,\n"
        "with site the site index of the origin of the figure."
      );

      bp::def
      (
        "find_pis",
        &find_pis2,
        ( bp::arg("classes"), bp::arg("structure"), bp::arg("site")=0 ),
        "Returns pis of a given structure for a given array of classes of equivalent figures,\n"
        "with site the site index of the origin of the figure."
      );
      bp::def
      (
        "find_pis2",
        &find_pis3,
        ( bp::arg("classes"), bp::arg("structure") ),
        "Returns pis of a given structure for a given array of classes of equivalent figures,\n"
        "with site the site index of the origin of the figure."
      );
    }

  }
} // namespace LaDa
