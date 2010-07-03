#include "LaDaConfig.h"

#include <boost/python/def.hpp>

#include "create_clusters.hpp"
#include "../create_clusters.h"


namespace LaDa
{
  namespace Python
  {
    boost::shared_ptr< CE::t_MLClusterClasses >
      create_clusters( Crystal::Lattice const &_lat, size_t _order, size_t _max, size_t _site )
        { return CE::create_clusters(_lat, _order, _max, _site ); }

    void expose_create_clusters()
    {
      namespace bp = boost::python;

      bp::def
      (
        "create_clusters",
        &create_clusters,
        (
          bp::arg("lattice"),
          bp::arg("order"),
          bp::arg("nth_shell"),
          bp::arg("site") = 0
        ),
        "Returns an array of arrays of equivalent clusters of given order, origin,"
        " and maximum distance (in number of shells) from the origin."
      );
    }

  }
} // namespace LaDa
