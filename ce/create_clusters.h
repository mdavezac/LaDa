//
//  Version: $Id$
//
#ifndef LADA_CE_CREATE_CLUSTERS_H
#define LADA_CE_CREATE_CLUSTERS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/shared_ptr.hpp>
#include "cluster.h"

namespace LaDa
{
  namespace CE
  {
    //! \brief Creates all clusters of oreder \a _order up to \a _n neighbors.
    void create_clusters( t_ClusterClasses &_out, Crystal::Lattice const &_lat,
                          size_t _order, size_t _neighbor, size_t _origin );

    //! \brief Creates all clusters of oreder \a _order up to \a _n neighbors.
    boost::shared_ptr<t_ClusterClasses> create_clusters( Crystal::Lattice const &_lat, 
                                                         size_t _order, size_t _neighbor, size_t _origin );

  } // end of namespace CE
} // namespace LaDa
#endif
