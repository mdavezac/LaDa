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
      create_clusters( Crystal::Lattice const &_lat, size_t _order, size_t _max, size_t _site )
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
        return CE::create_clusters(_lat, _order, _max, _site ); 
      }

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
