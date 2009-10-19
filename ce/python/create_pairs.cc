//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/shared_ptr.hpp>
#include <boost/python/def.hpp>
#include <boost/python/errors.hpp>

#include "create_pairs.hpp"
#include "../create_pairs.h"


namespace LaDa
{
  namespace Python
  {
    boost::shared_ptr< CE::t_ClusterClasses >
      create_pairs( Crystal::Lattice const &_lat, types::t_unsigned _max, size_t _site )
      {
        boost::shared_ptr< CE::t_ClusterClasses > result(new CE::t_ClusterClasses);
        CE::create_pairs(_lat, _max, *result, _site );
        return result;
      }

    void expose_create_pairs()
    {
      namespace bp = boost::python;

      bp::def
      (
       "create_pairs",
       &create_pairs,
       (
         bp::arg("lattice"),
         bp::arg("nbpairs"),
         bp::arg("site") = 0
       ),
       "Returns an array of nbpairs classes of equivalent pair clusters, with origin at _site in the lattice."
      );
    }

  }
} // namespace LaDa
