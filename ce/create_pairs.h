//
//  Version: $Id$
//
#ifndef LADA_CE_CREATE_PAIRS_H
#define LADA_CE_CREATE_PAIRS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>


#include "cluster.h"

namespace LaDa
{
  namespace CE
  {
    //! \brief Finds all pair clusters up to a given cutoff.
    //! \param[in]_lat is the input lattice.
    //! \param[in]_max_neigh is the number of pair classes.
    //! \param[inout] _out is the vector of pair classes.
    //! \param[in] index of the site for which to find pairs.
    void create_pairs( const Crystal :: Lattice &_lat,
                       types::t_unsigned _max_neigh,
                       std::vector< std::vector<Cluster> > &_out,
                       size_t _site = 0 );

  } // namespace CE

} // namespace LaDa

#endif
