//
//  Version: $Id$
//
#ifndef _CE_POSTOCONF_H_
#define _CE_POSTOCONF_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>
#include <opt/debug.h>

namespace Crystal
{
  class SplitIntoConfs
  {
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
      void operator()( const Crystal::Structure &_structure,
                       const types::t_unsigned _n );

      //! Returns a constant to the computed configurations.
      const t_Configurations& configurations() const { configurations_; }
      //! \brief Maps an integer position into an atomic index.
      size_t integer_to_index( size_t _atom ) const
      { 
        const size_t N( structure->atoms.size() );
        return o( ( ( _atom % (9*N) ) % (3*N) ) % N ); 
      }
      //! Clears results.
      void clear() { configurations_.clear(); structure = NULL; }


    protected: 
      //! Type of the container of position.
      typedef std::list< std::vector<size_t> > t_Positions;


      //! \brief Converts an integer position into an atomic positions.
      //! \details Natural numbers are used to map all atoms on the periodic
      //!          lattice and discriminate between periodic images. This
      //!          function goes from integer format to vectorial format.
      void integer_to_pos( size_t _atom, atat::rVector3d &_pos ) const;
      //! \brief Adds configurations, knowing the origin.
      //! \details Origin is indicated by first bit of \a _bitset::first on
      //!          input. The weight of the configuration will depend on \a
      //!          _bitset::second as given on input.  \a _bitset is modified
      //!          by this routine, except for the first bit of \a
      //!          _bitset::first and \a _bitset::second.
      void from_origin( t_CoefBitset &_bitset );
      //! \brief Compares positions indexed as integers, given the coordinate system.
      //! \note Does not compare norms of \a _a1 - \a _origin and \a _a2 - \a _origin,
      //!       E.g. does not implement complete comparison rule. 
      void compare_from_coords( const atat::rVector3d &_origin, 
                                const atat::rVector3d &_x, 
                                const atat::rVector3d &_y, 
                                const atat::rVector3d &_z, 
                                const atat::rVector3d &_a1, 
                                const atat::rVector3d &_a2 ) const; 
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
      const Crystal :: Structure * structure;
      //! The result.
      t_Configurations configurations_;
      //! The equivalent classes of positions.
      t_Positions epositions;
  };

} // end of Crystal namespace.

#endif
