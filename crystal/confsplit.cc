//
//  Version: $Id$
//

#include<algorithm>

#include<boost/bind.hpp>
#include<boost/ref.hpp>

#include<opt/algorithms.h>

#include "confsplit.h"

namespace Crystal
{
  void SplitIntoConfs :: operator()( const Crystal::Structure &_structure,
                                     const types::t_unsigned _n )
  {
    __DOASSERT( _n > 3, "Requested size of basis too small.\n" )
    n = _n;
    structure = &_structure;
    types::t_real weight( _structure.atoms.size() ); weight = 1e0 / weight;
    configurations_.clear();
    t_CoefBitset bitset( t_Bitset( _n ), weight );
    // This loop splits over possible origins.
    for( size_t i(0); i < _structure.atoms.size(); ++i )
    {
      bitset[0] = i;
      configurations_.push_back( bitset );
      from_origin();
    }
  }

  void integer_to_pos( size_t _atom, atat::rVector3d &_pos ) const
  {
    __ASSERT( not structure, "Structure pointer not set.\n" )
    const size_t N( structure->atoms.size() );
    const size_t i( _atom / (9*N) - 1 );
    const size_t j( ( _atom % (9*N) ) / (3*N) - 1 );
    const size_t k( ( ( _atom % (9*N) ) % (3*N) ) / N - 1 );
    const size_t o( ( ( _atom % (9*N) ) % (3*N) ) % N );
    _pos =   structure->cell *  atat::rVector3d( types::t_real(i),
                                                 types::t_real(j),
                                                 types::t_real(k) )
           + structure->atoms[o];
  }

  void SplitIntoConfs :: find_atoms_in_sphere( const atat::rVector3d &_origin,
                                               t_Positions &_positions )
  {
    __ASSERT( not structure, "Structure pointer not set.\n" )
    // first, computes and sorts nth neighbors.
    const size_t N( structure->atoms.size() );
    const size_t umax = N * ( ( (n / _structure->atoms.size() + 1) * 3 ) ** 3 + 1 );
    typedef std::pair< types::t_real, types::t_unsigned > t_normedpair;
    typedef std::vector< t_normedpair > t_normedpairs;
    t_normedpairs normedpairs; normedpairs.reseverve( umax );
    for( size_t u(0); u < umax; ++u )
    {
      atat::rVector3d dummy;
      integer_to_pos( u, dummy );
      dummy -= _origin;
      normedpair( std::pair< types::t_real, types::t_unsigned >( atat::norm2( dummy ), u ) );
    }
    std::partial_sort( normedpairs.begin(), normedpairs.begin() + n, normedpairs.end() );
    const types::t_real lastnorm( (normedpairs.begin() + n - 1)->first );


    // then  copies position as sets of eqivalent positions.
    // note that the std::partition makes sure we do not miss positions which
    // are the same distance as position n.
    _positions.clear();
    t_normedpairs::const_iterator i_pair = normedpairs.begin();
    t_normedpairs::const_iterator
      i_pair_end = std::partition
                   ( 
                     normedpairs.begin() + n, normedpairs.end(),
                     boost::bind( &t_normedpair::first, _1 ) < lastnorm
                   );
    lastnorm =  i_pair->first;
    while( i_pair != i_pair_end )
    {
      _positions.resize( _positions.size() + 1 );
      t_Positions :: value_type& back( _positions.back() );
      do
      {
        back.push_back( i_pair->second );
        ++i_pair 
      }
      while( i_pair != i_pair_end and Fuzzy( eq( i_pair->fist, lastnorm ) ) );
    }
  }

  void SplitIntoConfs :: from_origin( t_CoefBitset &_bitset )
  {
    __ASSERT( not structure, "Structure pointer not set.\n" )

    const atat::rVector3d origin( structure->atoms[ _bitset.first[0] ].pos );
    const types::t_real weight( _bitset.second );

    epositions.clear();
    find_atoms_in_sphere( origin, positions )
    // A copy of the positions used for last step sorting.
    t_Positions sorting_pos( positions );

    // loop over positions defining x.
    const t_Positions :: const_iterator i_xpositions = epositions.begin();
    const types::t_real xweight( _bitset.second / types::t_real( i_xpositions->size() ) );
    foreach( const t_Position::value_type::value_type &xindex, *i_xpositions )
    {
      atat::rVector3d dummy, dummy2; 
      integer( xindex, dummy );
      const atat::rVector3d x( dummy - origin );
      // finds positions defining y.
      const t_Positions :: const_iterator
         i_ypositions = ( i_xpositions->size() == 1 ) ? i_xpositions: i_xpositions + 1; 
      const types::t_real
        max_x_element = std::max_element
                        (
                           i_ypositions->begin(), i_ypositions->end(),
                           (
                             bl::bind( &SplitIntoConfs :: integer_to_pos,
                                       this, bl::_1, bl::var(dummy) ),
                             bl::bind( &SplitIntoConfs :: integer_to_pos,
                                       this, bl::_2, bl::var(dummy2) ),
                             bl::bind( &Fuzzy::le<types::t_real>,
                                       dummy * bl::constant(x), 
                                       dummy2 * bl::constant(x) )
                            )
                         );
      integer_to_pos( *maxelement, dummy ); 
      const types::t_real max_x_scalar_pos( dummy * x );
      // Stores possible y positions.
      std::vector< t_Positions::value_type::value_type > ypossibles;
      foreach( const t_Positions::value_type::value_type yindex, *i_ypositions )
      {
        integer_to_pos( yindex, dummy )
        if( Fuzzy::neq( dummy, max_x_scalar_pos ) ) continue;
        ypossibles.push_back( yindex );
      }

      _bitset.second = xweight / types::t_real( ypossibles.size() );
      foreach( const t_Positions::value_type::value_type yindex, ypossibles )
      {
        // at this point, we can define the complete coordinate system.
        integer( yindex, dummy );
        const atat::rVector3d y( dummy - origin );
        const atat::rVector3d z( x^y );

        // atoms are now included in the list according to the following rule:
        //  _ closest to the origin first.
        //  _ ties are broken according to largest x coordinate.
        //  _ next ties are broken according to largest y coordinate.
        //  _ final ties are broken according to largest z coordinate.

        // we iterate over distance from origin first.
        __ASSERT( _bitset.size() != n, "Incoherent sizes.\n" )
        t_CoefBitset::first::iterator i_bit = _bitset.first.begin();
        t_CoefBitset::first::iterator i_bit_end = _bitset.first.end();
        size_t nbit(0);
        foreach( t_Positions :: value_type equaldistance, sorting_pos )
        {
          if( nbit == n ) break;
          __ASSERT( nbit > n, "index out of range.\n" )

          const size_t edn( std::min( equaldistance.size(), n - nbit ) );
          if( edn == 1 ) 
          {
            __ASSERT( nbit == 0 and not ( *i_bit != equaldistance.front() ),
                      "Incoherent sorted position and origin.\n" )
            __ASSERT( nbit == 1 and not ( yindex != equaldistance.front() ),
                      "Incoherent sorted position and y coord.\n" )
            *i_bit = equaldistance.front();
            ++i_bit;
            ++nbit;
            continue;
          }

          if( edn <= n - nbit ) 
            std::sort
            ( 
              sorting_pos.begin(), sorting_pos.end(),
              boost::bind( &SplitIntoConfs::compare_from_coords,
                           this, origin, x, y, z, _1, _2 );
            );
          else  std::partial_sort
                ( 
                  sorting_pos.begin(), sorting_pos.begin() + n - nbit, sorting_pos.end(),
                  boost::bind( &SplitIntoConfs::compare_from_coords, this, _1, x, y, z );
                );
          __ASSERT( nbit == 0 and not ( *i_bit != equaldistance.front() ), 
                      "Incoherent sorted position and origin.\n" )
          __ASSERT(     nbit <= 1 and  nbit + edn > 1 
                    and not ( yindex != equaldistance[ 1 - nbit ] ), 
                      "Incoherent sorted position and y coord.\n" )
          std::copy( equaldistance.begin(), equaldistance.begin() + edn,
                     boost::ref( i_bit ) );
          nbit += edn;
        } // end of loop over positions at equivalent distance.


        // finally adds configuration.
        configurations_.push_back( bitset );
      } // end of loop over equivalent y coords.
  
    } // end of loop over equivalent  x coords.

    // resets weight.
    _bitset.second = weight;
  }

    void compare_from_coords( const atat::rVector3d &_origin, 
                              const atat::rVector3d &_x, 
                              const atat::rVector3d &_y, 
                              const atat::rVector3d &_z, 
                              const atat::rVector3d &_a1, 
                              const atat::rVector3d &_a2 ) const
    { 
      const atat::rVector3d b1( _a1 - _origin );
      const atat::rVector3d b2( _a2 - _origin );

      __ASSERT( Fuzzy::neq( atat::norm2( b1 ), atat::norm2( b2 ) ),
                "Expected distance to origin to be equal.\n" )

      const types::t_real x1( b1 * _x ), x2( b2 * _x );
      if( Fuzzy::neq( x1, x2 ) ) return Fuzzy::gt( x1, x2 );
      const types::t_real y1( b1 * _y ), y2( b2 * _y );
      if( Fuzzy::neq( y1, y2 ) ) return Fuzzy::gt( y1, y2 );
      return Fuzzy::gt( b1 * _z, b2 * _z );
    }
  }
} // namespace Crystal

