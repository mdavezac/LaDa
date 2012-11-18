#ifndef LADA_CRYSTAL_MAPSITES_H
#define LADA_CRYSTAL_MAPSITES_H

#include "LaDaConfig.h"

#include "structure/structure.h"

namespace LaDa 
{
  namespace crystal
  {
    namespace
    {
      //! \brief Map atomic sites from mapper onto mappee.
      //! \param[in] _mapper : a lattice against which to map atomic sites.
      //! \param[inout] _mappee : a supercell for which sites will be mapped.
      //! \param[in] _withocc : whether to take occupation into count, or only position.
      //! \param[in] _tolerance : Tolerance criteria for distances, in units of _mappee.scale().
      //! \details Where possible, the site indices of the mappee structure
      //!          corresponds to the equivalent sites in the mapper structure.
      //! \return True if mapping is successful, False if all sites could not be mapped. 
      //!         Since in the case of defects, incomplete mappings may be what is wanted, 
      static bool map_sites( Structure const &_mapper, Structure &_mappee,
                             python::Object _withocc = python::Object(),
                             types::t_real _tolerance = types::tolerance )
#     ifdef LADA_CRYSTAL_MODULE
          ;
#     else
          { return (*(bool(*)(Structure const&, Structure &, python::Object, types::t_real))
                    api_capsule[10])(_mapper, _mappee, _withocc, _tolerance); }
#     endif
    }
  } // namespace crystal
} // namespace LaDa

#endif
