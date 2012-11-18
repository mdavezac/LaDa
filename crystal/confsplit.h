#ifndef _CRYSTAL_CONFSPLIT_H_
#define _CRYSTAL_CONFSPLIT_H_

#include "LaDaConfig.h"

#include <vector>
#include <utility>

#include <misc/types.h>

#include "structure/structure.h"

namespace LaDa
{
  namespace crystal
  {
    namespace // limits declaration to current translation unit.
    {
      //! \brief Creates a split-configuration for a given structure and atomic origin.
      //! \details Split-configurations are a symmetry-agnostic atom-centered
      //!          description of chemical environment. For details, see
      //!          \url{http://dx.doi.org/10.1137/100805959} "d'Avezac, Botts,
      //!          Mohlenkamp, Zunger, SIAM J. Comput. \Bold{30} (2011)". 
      //! \param[in] _structure Structure for which to create split-configurations.
      //! \param[in] _index Index into _structure for the origin of the configuration.
      //! \param[in] _nmax Number of atoms (cutoff) to consider for inclusion
      //!                  in the split-configuration.
      //! \param[inout] _configurations Object where configurations should be
      //!       stored. It should a null object,  or a list of previously
      //!       existing configurations. There is no error checking, so do not
      //!       mix and match.
      //! \param[in] _tolerance Tolerance criteria when comparing distances.
      //! \return A list of splitted configuration. Each item in this list is
      //!         itself a list with two inner items. The first inner item is an
      //!         ordered list of references to atoms. The second inner item is
      //!         the weight for that configuration. The references to the atoms
      //!         are each a 3-tuple consisting of an actual reference to an
      //!         atom, a translation vector from the center of the configuration
      //!         to the atom's relevant periodic image, and a distance from the
      //!         center. [[[(atom, vector from center, distance from center),
      //!         ...], weight], ...]
      bool splitconfigs( Structure const &_structure,
                         Atom const &_origin,
                         Py_ssize_t _nmax,
                         python::Object &_configurations,
                         types::t_real _tolerance )
#     ifdef LADA_CRYSTAL_MODULE
        ;
#     else
          { return (*(bool(*)( Structure const&, Atom const&, 
                               Py_ssize_t, python::Object &,
                               types::t_real))
                    api_capsule[18])(_structure, _origin, _nmax, _configurations, _tolerance); }
#     endif
    } // anonymous namespace.
  } // end of Crystal namespace.
} // namespace LaDa

#endif
