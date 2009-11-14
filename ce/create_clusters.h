//
//  Version: $Id$
//
#ifndef LADA_CE_CREATE_CLUSTERS_H
#define LADA_CE_CREATE_CLUSTERS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/shared_ptr.hpp>
#include "mlclusters.h"

namespace LaDa
{
  namespace CE
  {
    //! \brief Creates all clusters of order \a _order up to \a _n neighbors.
    void create_clusters( t_MLClusterClasses &_out, Crystal::Lattice const &_lat,
                          size_t const _order, size_t const _neighbor, size_t const _origin = 0 );
    //! \brief Creates all clusters of order \a _order up to \a _n neighbors.
    boost::shared_ptr<t_MLClusterClasses>
      create_clusters(Crystal::Lattice const &_lat, size_t const _order,
                      size_t const _neighbor, size_t const _origin = 0);
  } // end of namespace CE
} // namespace LaDa
#endif
