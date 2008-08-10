//
//  Version: $Id: separable.h 626 2008-06-18 22:05:05Z davezac $
//
#ifndef _CE_POSTOCONF_H_
#define _CE_POSTOCONF_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <vector>

#include <opt/types.h>
#include <atat/vectmac.h>
#include <crystal/lattice.h>
#include <crystal/structure.h>

#include "separable.h"
#include "boolean.h"

namespace CE
{
  //! Prepares a fixed-lattice separable function from input.
  class PosToConfs
  {
    public:
      //! Type of the container of positions.
      typedef std::vector< atat::rVector3d > t_Basis;
      //! Type of the container of symmetry operations.
      typedef std::vector< atat::rMatrix3d > t_SymOps;
      //! Type of the pure bitset representing a configuration for a set symmetry.
      typedef std::vector< bool > t_Bitset;
      //! Type  containing a pure bitset and an attached coefficient.
      typedef std::pair< t_Bitset, types::t_real > t_CoefBitset;
      //! Type representing a set of configurations with all possible symmetries.
      typedef std::vector< t_CoefBitset > t_Configurations;

      //! the positions as currently computed.
      t_Basis basis;

      //! Constructor.
      PosToConfs ( Crystal::Lattice &_lat ) { init_syms( _lat ); }
      //! Creates a basis of positions.
      void create_positions( const std::string &_desc );
      
      //! Creates all necessary configurations for a given structure.
      void operator()( const Crystal :: Structure &_structure,
                       t_Configurations& _confs ) const;

    protected:
      //! Initializes list of symmetry operations.
      void init_syms ( Crystal::Lattice &_lat );
      
      //! The syemmetry operations.
      t_SymOps syms;
  };

} // end of CE namespace

#endif //  _SEPARABLE_CE_H_
