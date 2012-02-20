#include "LaDaConfig.h"

#include <boost/python/def.hpp>
#include <boost/python/errors.hpp>

#include "create_clusters.hpp"
#include "../create_clusters.h"


namespace LaDa
{
  namespace Python
  {
    boost::shared_ptr< CE::t_MLClusterClasses >
      create_clusters( Crystal::Lattice &_lat, size_t _order, size_t _max, size_t _site )
      {
        Crystal::Lattice :: t_Sites :: const_iterator i_first = _lat.sites.begin();
        Crystal::Lattice :: t_Sites :: const_iterator const i_end = _lat.sites.end();
        bool found = false;
        for(; i_first != i_end and not found; ++i_first) found = i_first->type.size() > 1;
        if(not found) 
        {
          PyErr_SetString( PyExc_ValueError,
                           "Lattice never has more thant one specie per site. "
                           "Cannot perform CE." );
          boost::python::throw_error_already_set();
          return boost::shared_ptr<CE::t_MLClusterClasses>();
        }
        Crystal::Lattice *old0 = Crystal::Structure::lattice;
        Crystal::Lattice *old1 = Crystal::TStructure<std::string>::lattice;
        Crystal::Structure::lattice = &_lat;
        Crystal::TStructure<std::string>::lattice = &_lat;
        boost::shared_ptr<CE::t_MLClusterClasses> result( CE::create_clusters(_lat, _order, _max, _site ) ); 
        Crystal::Structure::lattice = old0;
        Crystal::TStructure<std::string>::lattice = old1;
        return result;
      }

    void expose_create_clusters()
    {
      namespace bp = boost::python;

      bp::def
      (
        "_create_clusters",
        &create_clusters,
        (
          bp::arg("lattice"),
          bp::arg("order"),
          bp::arg("nth_shell"),
          bp::arg("site") = 0
        ),
        "Creates arrays of equivalence clusters for given order, origin, and distance.\n\n"
        ":Parameters:\n"
        "  lattice\n    Back-bone lattice of the Ising-model.\n" 
        "  order : int\n    Size of the cluster of interest, eg "\
            "on-site, pair, or higher-order interaction.\n"
        "  nth_shell\n    Clusters include spins up to nth shell removed (inclusive).\n" 
        "  site\n    Origin of the cluster in the case of multi-site lattices.\n" 
      );
    }

  }
} // namespace LaDa
