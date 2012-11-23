#if LADA_CRYSTAL_MODULE != 1

#include "LaDaConfig.h"

namespace LaDa 
{
  namespace crystal
  {
    namespace
    {
#endif

#if LADA_CRYSTAL_MODULE != 1
  //! \brief Map atomic sites from mapper onto mappee.
  //! \param[in] _mapper : a lattice against which to map atomic sites.
  //! \param[inout] _mappee : a supercell for which sites will be mapped.
  //! \param[in] _withocc : whether to take occupation into count, or only position.
  //! \param[in] _tolerance : Tolerance criteria for distances, in units of _mappee.scale().
  //! \details Where possible, the site indices of the mappee structure
  //!          corresponds to the equivalent sites in the mapper structure.
  //! \return True if mapping is successful, False if all sites could not be mapped. 
  //!         Since in the case of defects, incomplete mappings may be what is wanted, 
  bool map_sites( Structure const &_mapper, Structure &_mappee,
                  python::Object _withocc = python::Object(),
                  types::t_real _tolerance = types::tolerance )
    LADA_END( { return (*(bool(*)(Structure const&, Structure &, python::Object, types::t_real))
                        api_capsule[LADA_SLOT(crystal)])(_mapper, _mappee, _withocc, _tolerance); } )
#else
  api_capsule[LADA_SLOT(crystal)] = (void *)map_sites;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(crystal))
#include LADA_ASSIGN_SLOT(crystal)

#if LADA_CRYSTAL_MODULE != 1
    }
  } // namespace crystal
} // namespace LaDa
#endif
