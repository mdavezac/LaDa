//
//  Version: $Id$
//
#ifndef _CE_POSTOCONF_H_
#define _CE_POSTOCONF_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <vector>
#include <utility>

#include <opt/types.h>
#include <atat/vectmac.h>
#include <crystal/lattice.h>
#include <crystal/structure.h>

namespace LaDa
{

  namespace CE
  {
    //! Prepares a fixed-lattice separable function from input.
    class PosToConfs
    {
      public:
        //! Type of the container of positions.
        typedef std::vector< atat::rVector3d > t_Positions;
        //! Type of the container of symmetry operations.
        typedef std::vector< atat::rMatrix3d > t_SymOps;
        //! Type of the pure bitset representing a configuration for a set symmetry.
        typedef std::vector< types::t_unsigned > t_Bitset;
        //! Type  containing a pure bitset and an attached coefficient.
        typedef std::pair< t_Bitset, types::t_real > t_CoefBitset;
        //! Type representing a set of configurations with all possible symmetries.
        typedef std::vector< t_CoefBitset > t_Configurations;

        //! the positions as currently computed.
        t_Positions positions;

        //! Constructor.
        PosToConfs ( Crystal::Lattice &_lat = *Crystal::Structure::lattice ) : n(0, 0)
          { init_syms( _lat ); }
        //! Copy Constructor.
        PosToConfs ( const PosToConfs &_c ) : syms( _c.syms ), n(_c.n) {}
        //! Creates a basis of positions.
        void create_positions( const std::string &_desc );
        
        //! Creates all necessary configurations for a given structure.
        void operator()( const Crystal :: Structure &_structure,
                         t_Configurations& _confs ) const;

        //! Creates all necessary configurations for a given structure for position basis.
        void posbasis( const Crystal :: Structure &_structure,
                       t_Configurations& _confs ) const;

        //! Creates all necessary configurations for a given structure for position n sites.
        void nsites( const Crystal :: Structure &_structure,
                     t_Configurations& _confs ) const;
        //! Number of degrees of freedom.
        size_t dof() const
          { return n.second - n.first ? n.second - n.first+1: positions.size(); }

      protected:
        //! Initializes list of symmetry operations.
        void init_syms ( Crystal::Lattice &_lat );
        
        //! The symmetry operations.
        t_SymOps syms;
        //! number of sites.
        std::pair<size_t, size_t> n;
    };

    inline void PosToConfs :: operator()( const Crystal :: Structure &_structure,
                                          t_Configurations& _confs ) const
    {
      __ASSERT( n.second - n.first == 0 and positions.empty(), "Nothing initialized.\n" )
      ( n.second - n.first == 0 ) ?
        posbasis( _structure, _confs ): 
        nsites( _structure, _confs );
    }

  } // end of CE namespace
} // namespace LaDa
#endif //  _SEPARABLE_CE_H_
