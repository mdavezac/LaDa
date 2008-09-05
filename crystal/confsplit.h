//
//  Version: $Id$
//
#ifndef _CRYSTAL_CONFSPLIT_H_
#define _CRYSTAL_CONFSPLIT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <utility>

#include <opt/types.h>
#include <opt/debug.h>

#include "structure.h"

namespace Crystal
{
  class SplitIntoConfs
  {
    protected: 
      //! Type of the of a position.
      typedef std::pair< atat::rVector3d, size_t > t_Position;
      //! Type of the container of positions.
      typedef std::list< std::vector< t_Position > > t_Positions;

    public:
      //! Type of the pure bitset representing a configuration for a set symmetry.
      typedef std::vector< types::t_unsigned > t_Bitset;
      //! Type  containing a pure bitset and an attached coefficient.
      typedef std::pair< t_Bitset, types::t_real > t_CoefBitset;
      //! Type representing a set of configurations with all possible symmetries.
      typedef std::vector< t_CoefBitset > t_Configurations;

      //! Constructor
      SplitIntoConfs() : n(0), structure(NULL) {}

      //! Computes the configurations.
      void operator()( const Structure &_structure,
                       const types::t_unsigned _n );

      //! Returns a constant to the computed configurations.
      const t_Configurations& configurations() const { configurations_; }
      //! Clears results.
      void clear() { configurations_.clear(); structure = NULL; }


    protected: 
      //! \brief Adds configurations, knowing the origin.
      //! \details Origin is indicated by first bit of \a _bitset::first on
      //!          input. The weight of the configuration will depend on \a
      //!          _bitset::second as given on input.  \a _bitset is modified
      //!          by this routine, except for the first bit of \a
      //!          _bitset::first and \a _bitset::second.
      void from_origin( t_CoefBitset &_bitset );
      //! \brief Compares positions indexed as integers, given the origin and x-axis.
      bool compare_from_x( const atat::rVector3d &_origin, 
                           const atat::rVector3d &_x, 
                           const t_Position& _a1, 
                           const t_Position& _a2 ) const
        { return Fuzzy::le( _a1.first * _x, _a2.first * _x ); }
      //! \brief Compares positions indexed as integers, given the coordinate system.
      //! \note Does not compare norms of \a _a1 - \a _origin and \a _a2 - \a _origin,
      //!       E.g. does not implement complete comparison rule. 
      bool compare_from_coords( const atat::rVector3d &_origin, 
                                const atat::rVector3d &_x, 
                                const atat::rVector3d &_y, 
                                const atat::rVector3d &_z, 
                                const t_Position& _a1, 
                                const t_Position& _a2 ) const; 
      //! \brief Computes a list of atoms of the SplitIntoCoefs::n closest
      //!        neighbors to the origin.
      //! \params[in] _origin from which to compute distances.
      //! \params[out] _positions a vector of vectors positions which are at
      //!                            the same distance from the origin. E.g.
      //!                            each internal vector contains all
      //!                            positions on  a sphere.
      void find_atoms_in_sphere( const atat::rVector3d &_origin,
                                 t_Positions &_positions );

      //! The number of atoms to included in each configuration.
      types::t_unsigned n;
      //! Structure for which configurations are computed.
      const  Structure * structure;
      //! The result.
      t_Configurations configurations_;
      //! The equivalent classes of positions.
      t_Positions epositions;
  };

} // end of Crystal namespace.

#endif
