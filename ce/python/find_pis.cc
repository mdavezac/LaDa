#include "LaDaConfig.h"

#include <iostream>
#include <sstream>
#include <complex>

#include <boost/python/def.hpp>
#include <boost/python/list.hpp>


#include <python/numpy_types.h>
#include "find_pis.hpp"
#include "../find_pis.h"


namespace LaDa
{
  namespace Python
  {
    namespace bp = boost::python;
    bp::object find_pis( CE::t_ClusterClasses const &_cls,
                         Crystal::TStructure<std::string> const &_str,
                         size_t _site )
    {
      Crystal::Structure str;
      Crystal::convert_string_to_real_structure(_str, str);
      std::vector<types::t_real> vec;
      CE::find_pis(_cls, str, vec, _site);
      return math::numpy::copy_1darray(vec);
    }
    bp::object find_pis3( CE::t_ClusterClasses const &_cls,
                          Crystal::TStructure<std::string> const &_str)
    {
      Crystal::Structure str;
      Crystal::convert_string_to_real_structure(_str, str);
      std::vector<types::t_real> vec;
      CE::find_pis2(_cls, str, vec);
      return math::numpy::copy_1darray(vec);
    }

//   bp::object find_pis2( bp::list const &_cls,
//                         Crystal::TStructure<std::string> const &_str,
//                         size_t _site )
//   {
//     CE::t_ClusterClasses cls(bp::len(_cls));
//     for(size_t i(0); i < bp::len(_cls); ++i)
//       cls[i] = bp::extract<CE::t_Clusters>(_cls[i]);
//     return find_pis(cls, _str, _site);
//   }

    void expose_find_pis()
    {
      import_array(); // needed for NumPy 

      bp::def( "find_pis", &find_pis, 
               (bp::arg("classes"), bp::arg("structure"), bp::arg("site")=0) );
 //    bp::def
 //    (
 //      "find_pis",
 //      &find_pis2,
 //      ( bp::arg("classes"), bp::arg("structure"), bp::arg("site")=0 ),
 //      "Computes S{Pi}s.\n\n"
 //      "@param classes: classes of quivalents cluster classes.\n"
 //      "@param structure: structure for which to compute S{Pi}s.\n"
 //      "@type structure: L{crystal.Structure}\n"
 //      "@param site: Which site to compute S{Pi}s for. Default 0.\n"
 //      "@return: S{Pi}s for given structure and cluster classes.\n"
 //    );
      bp::def
      (
        "find_pis2",
        &find_pis3,
        ( bp::arg("classes"), bp::arg("structure") ),
        "Computes S{Pi}s.\n\n"
        "@param classes: classes of quivalents cluster classes.\n"
        "@param structure: structure for which to compute S{Pi}s.\n"
        "@type structure: L{crystal.Structure}\n"
        "@param site: Which site to compute S{Pi}s for. Default 0.\n"
        "@return: S{Pi}s for given structure and cluster classes.\n"
      );
    }

  }
} // namespace LaDa
