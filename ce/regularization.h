//
//  Version: $Id: functional_builder.impl.h 685 2008-07-22 17:20:39Z davezac $
//

#ifndef _CE_REGULARIZATION_H_
#define _CE_REGULARIZATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>
#include <opt/debug.h>
#include <crystal/structures.h>

namespace CE
{
  //! \brief Computes pis of \a _str for \a _clusters.
  //! \param[in] _cluster is a vector of containers of symmetrically equivalent
  //!                     cluster, centered on the origin.
  //! \param[in] _str the structure for which to compute the pis.
  //! \param[out] _pis the computed pis, one per class of symmetrically
  //!                  equivalent clusters.
  void find_pis( const std::vector< std::vector< Cluster > > &_clusters,
                 const Crystal::Structure & _str,
                 std::vector< types::t_int > &_pis );
  //! \brief Computes pis of \a _str for \a _clusters.
  //! \see[in] _cluster is a vector of containers of symmetrically equivalent
  //!                     cluster, centered on the origin.
  //! \param[in] _str structures for which to compute the pis.
  //! \param[out] _pis the computed pis, one per class of symmetrically
  //!                  equivalent clusters.
  void find_pis( const std::vector< std::vector< Cluster > > &_clusters,
                 std::vector< Crystal::Structure > & _str,
                 std::vector< std::vector< types::t_int > > &_pis );
} // end of namespace CE

#endif 
