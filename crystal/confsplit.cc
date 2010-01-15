//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<algorithm>

#include<boost/bind.hpp>
#include<boost/ref.hpp>
#include<boost/lambda/bind.hpp>
// #include<boost/math/special_functions/pow.hpp>

#include<opt/algorithms.h>

#include "confsplit.h"

namespace LaDa
{
  namespace Crystal
  {
    //! \cond
    namespace details
    {
      template<typename _ForwardIterator, typename _Compare, typename _Skip >
      _ForwardIterator
      max_element(_ForwardIterator __first, _ForwardIterator __last,
                  _Compare __comp, _Skip __skip )
      {
        // concept requirements
    //   __glibcxx_function_requires(_ForwardIteratorConcept<_ForwardIterator>)
    //   __glibcxx_function_requires(_BinaryPredicateConcept<_Compare,
    //         typename iterator_traits<_ForwardIterator>::value_type,
    //         typename iterator_traits<_ForwardIterator>::value_type>)
    //   __glibcxx_requires_valid_range(__first, __last);

        if (__first == __last) return __first;
        do if( not __skip( *__first ) ) break;
        while (++__first != __last);
        if( __first == __last ) return __last;
        _ForwardIterator __result = __first;
        while (++__first != __last)
        {
          if (__skip( *__first ) ) continue;
          if (__comp(*__result, *__first)) __result = __first;
        }
        return __result;
      }
    }
    //! \endcond

    void SplitIntoConfs :: operator()( const Structure &_structure,
                                       const size_t _n )
    {
      __DOASSERT( _n <= 3, "Requested size of basis too small.\n" )
      n.first = 0; n.second = _n;
      structure = &_structure;
      configurations_.clear();
      // This loop splits over possible origins.
      for( size_t i(0); i < _structure.atoms.size(); ++i )
        from_origin( i );
    }
    void SplitIntoConfs :: operator()( const Structure &_structure,
                                       const std::pair< size_t, size_t > _n )
    {
      __DOASSERT( _n.second - _n.first <= 3, "Requested size of basis too small.\n" )
      n = _n;
      structure = &_structure;
      configurations_.clear();
      // This loop splits over possible origins.
      for( size_t i(0); i < _structure.atoms.size(); ++i )
        from_origin( i );
    }

    void SplitIntoConfs :: find_atoms_in_sphere( const Eigen::Vector3d &_origin,
                                                 t_Positions &_positions )
    {
      //namespace bm = boost::math;
      __ASSERT( not structure, "Structure pointer not set.\n" )
      // first, computes and sorts nth neighbors.
      const types::t_int N( structure->atoms.size() );
      const types::t_int umax = n.second / structure->atoms.size() + 1;
      typedef std::pair< types::t_real, t_Position > t_normedpair;
      typedef std::vector< t_normedpair > t_normedpairs;
      t_normedpairs normedpairs; normedpairs.reserve( N * ( umax * 4 )^3 );
      for( types::t_int x(-umax); x <= umax; ++x )
        for( types::t_int y(-umax); y <= umax; ++y )
          for( types::t_int z(-umax); z <= umax; ++z )
          {
            Structure :: t_Atoms :: const_iterator i_atom = structure->atoms.begin();
            Structure :: t_Atoms :: const_iterator i_atom_end = structure->atoms.end();
            t_normedpair pair;
            const Eigen::Vector3d trans(   structure->cell
                                         * Eigen::Vector3d( types::t_real(x),
                                                            types::t_real(y),
                                                            types::t_real(z) ) );
            for( pair.second.second=0; i_atom != i_atom_end;
                 ++i_atom, ++pair.second.second )
            {
              pair.second.first = trans + i_atom->pos - _origin; 
              pair.first = pair.second.first.squaredNorm();
              normedpairs.push_back( pair );
            }
          }
      std::partial_sort
      ( 
        normedpairs.begin(), normedpairs.begin() + n.second, normedpairs.end(),
        boost::bind
        (
          &Fuzzy::le<types::t_real>, 
          boost::bind( &t_normedpair::first, _1), 
          boost::bind( &t_normedpair::first, _2) 
        )
      );


      // then  copies position as sets of equivalent positions.
      // note that the std::partition makes sure we do not miss positions which
      // are the same distance as position n.
      _positions.clear();
      const types::t_real lastnorm( (normedpairs.begin() + n.second - 1)->first );
      t_normedpairs::const_iterator i_pair = normedpairs.begin();
      t_normedpairs::const_iterator
        i_pair_end = std::partition
                     ( 
                       normedpairs.begin() + n.second, normedpairs.end(),
                       boost::bind
                       ( 
                         &Fuzzy::eq<types::t_real>,
                         boost::bind( &t_normedpair::first, _1 ),
                         lastnorm
                       )
                     );
      size_t i( n.first );
      for(; i_pair != i_pair_end and i > 0; ++i_pair )
      {
        if( Fuzzy::is_zero(i_pair->first) ) continue; 
        --i;
      }
      __DOASSERT( i_pair == i_pair_end, "Insufficient number of atoms considered.\n" )

      while( i_pair != i_pair_end )
      {
        if( Fuzzy::is_zero(i_pair->first) ) { ++i_pair; continue; }

        const types::t_real current_norm = i_pair->first;
        _positions.resize( _positions.size() + 1 );
        t_Positions :: value_type& back( _positions.back() );
        do
        {
          back.push_back( i_pair->second );
          ++i_pair;
        }
        while( i_pair != i_pair_end and Fuzzy::eq( i_pair->first, current_norm ) );
      }
    }

    void SplitIntoConfs :: from_origin( size_t _i )
    {
      namespace bl = boost::lambda;
      __ASSERT( not structure, "Structure pointer not set.\n" )

      const types::t_real weight( 1e0 / types::t_real(structure->atoms.size()) );
      const size_t nsize( n.second - n.first );
      t_CoefBitset bitset( t_Bitset( nsize  ), weight );
      bitset.first[0] = t_Position( structure->atoms[_i].pos, _i );
      const Eigen::Vector3d origin( structure->atoms[ _i ].pos );

      epositions.clear();
      find_atoms_in_sphere( origin, epositions );
      // A copy of the epositions used for last step sorting.
      t_Positions sorting_pos( epositions );

      // loop over epositions defining x.
      const t_Positions :: const_iterator i_xpositions = epositions.begin();
      const types::t_real
        xweight( bitset.second / types::t_real( i_xpositions->size() ) );
      foreach( const t_Positions::value_type::value_type &xPos, *i_xpositions )
      {
        const Eigen::Vector3d x( xPos.first - origin );
        // finds positions defining y.
        // Stores possible y positions.
        std::vector< t_Position > ypossibles;
        t_Positions :: const_iterator i_ypositions = i_xpositions;
        if( i_xpositions->size() == 1 ) ++i_ypositions; 
        // Pointers are defined explicitely as a workaround for commercial
        // compilers, such as pgi.
        bool (*ptr_is_zero)( types::t_real ) = &Fuzzy::is_zero<types::t_real>;
        const t_Positions :: value_type :: const_iterator 
          max_x_element = details::max_element
                          (
                            i_ypositions->begin(), i_ypositions->end(),
                            boost::bind( &SplitIntoConfs::compare_from_x,
                                         this, origin, x, _1, _2 ),
                            bl::bind
                            (
                              ptr_is_zero,
                              bl::bind
                              ( 
                                &Eigen::Vector3d::squaredNorm,
                                  bl::bind( &t_Position::first, bl::_1 ) 
                                - bl::constant(xPos.first)
                              )
                            )
                          );
        if( max_x_element == i_ypositions->end() ) 
        {
          std::cout << *structure << "\n";
          foreach( const t_Position &_pos, *i_ypositions )
            std::cout << _pos.first << "\n";
          __DOASSERT( true, i_ypositions->size() << "\n" )
        }
        const types::t_real max_x_scalar_pos( max_x_element->first * x );
        foreach( const t_Position yPos, *i_ypositions )
        {
          if( Fuzzy::neq( yPos.first * x, max_x_scalar_pos ) ) continue;
          if( Fuzzy::is_zero( (yPos.first - xPos.first).squaredNorm() ) ) continue;
          ypossibles.push_back( yPos );
        }

        bitset.second = xweight / types::t_real( ypossibles.size() );
        foreach( const t_Position yPos, ypossibles )
        {
          // at this point, we can define the complete coordinate system.
          const Eigen::Vector3d y( yPos.first - origin );
          const Eigen::Vector3d z( x^y );

          // atoms are now included in the list according to the following rule:
          //  _ closest to the origin first.
          //  _ ties are broken according to largest x coordinate.
          //  _ next ties are broken according to largest y coordinate.
          //  _ final ties are broken according to largest z coordinate.

          // we iterate over distance from origin first.
          __ASSERT( bitset.first.size() != nsize, "Bitset too small.\n" );
          t_CoefBitset::first_type::iterator i_bit = bitset.first.begin();
          t_CoefBitset::first_type::iterator i_bit_end = bitset.first.end();
          size_t nbit(1);
          foreach( t_Positions :: value_type equaldistance, sorting_pos )
          {
            if( nbit == nsize ) break;
            __ASSERT( nbit > nsize, "index out of range.\n" )

            const size_t edn( std::min( equaldistance.size(), nsize - nbit ) );
            if( edn == 1 ) 
            {
              *i_bit = equaldistance.front();
              ++i_bit;
              ++nbit;
              continue;
            }

            if( edn <= nsize - nbit ) 
              std::sort
              ( 
                equaldistance.begin(), equaldistance.end(),
                boost::bind( &SplitIntoConfs::compare_from_coords,
                             this, origin, x, y, z, _1, _2 )
              );
            else std::partial_sort
                 ( 
                   equaldistance.begin(),
                   equaldistance.begin() + nsize - nbit, 
                   equaldistance.end(),
                   boost::bind( &SplitIntoConfs::compare_from_coords,
                                this, origin, x, y, z, _1, _2 )
                 );
            std::copy( equaldistance.begin(), equaldistance.begin() + edn, i_bit );
            i_bit += edn;
            nbit += edn;
          } // end of loop over positions at equivalent distance.


          // finally adds configuration.
          t_Configurations :: iterator i_found = configurations_.begin();
          t_Configurations :: iterator i_conf_end = configurations_.end();
          for(; i_found != i_conf_end; ++i_found )
            if( bitset.first == i_found->first ) break;
          if( i_found == i_conf_end ) configurations_.push_back( bitset );
          else i_found->second += bitset.second;
        } // end of loop over equivalent y coords.
      } // end of loop over equivalent  x coords.
    }

    bool SplitIntoConfs :: compare_from_coords( const Eigen::Vector3d &_origin, 
                                                const Eigen::Vector3d &_x, 
                                                const Eigen::Vector3d &_y, 
                                                const Eigen::Vector3d &_z, 
                                                const t_Position& _a1, 
                                                const t_Position& _a2 ) const
    { 
      const types::t_real x1( _a1.first * _x ), x2( _a2.first * _x );
      if( Fuzzy::neq( x1, x2 ) ) return Fuzzy::gt( x1, x2 );
      const types::t_real y1( _a1.first * _y ), y2( _a2.first * _y );
      if( Fuzzy::neq( y1, y2 ) ) return Fuzzy::gt( y1, y2 );
      return Fuzzy::gt( _a1.first * _z, _a2.first * _z );
    }



  } // namespace Crystal
} // namespace LaDa
