//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <algorithm>
#include <functional>

#include <opt/ndim_iterator.h>
#include <opt/atat.h>

#include "vff.h"
  
namespace LaDa
{
  namespace Vff
  { 
    namespace details
    {
      template<  class T_CONTAINER >
        class Order
        {
          public:
            bool operator( atat::rVector3d _a, atat::rVector3d _b ) const
            {
              const types::t_real a( atat::norm2(_a) );
              const types::t_real b( atat::norm2(_b) );
              if( Fuzzy::neq(a, b) ) return a < b;
              if( Fuzzy::neq( _a(0), _b(0) )  ) return _a(0) < _b(0);
              if( Fuzzy::neq( _a(1), _b(1) )  ) return _a(1) < _b(1);
              return _a(2) < _b(2);
            }
        };
      template<  class T_CONTAINER >
        class OrderI
        {
          public:
            OrderI  ( atat::rVector3d &_origin, 
                      const T_CONTAINER &_container )
                   : o( _origin ), c( _container ) {}
            OrderI( const Order &_c ) : o(_c.o), c(_c.c) {}

            bool operator( atat::rVector3d _a, atat::rVector3d _b ) const
            {
              const types::t_real a( atat::norm2(_a) );
              const types::t_real b( atat::norm2(_b) );
              if( Fuzzy::neq(a, b) ) return a < b;
              if( Fuzzy::neq( _a(0), _b(0) )  ) return _a(0) < _b(0);
              if( Fuzzy::neq( _a(1), _b(1) )  ) return _a(1) < _b(1);
              return _a(2) < _b(2);
            }
            bool operator( size_t _a, size_t _b ) const;
            {
              ASSERT( _a < c.size() and _b < c.size(), "Index out-of-range.\n" )
              return Order(  c[_a] - o,  c[_b] - o );
            }
          private:;
            const T_CONTAINER &c;
            const atat::rVector3d &o;
        };
    }

    bool Vff :: build_tree_sort()
    {
      // finds first neighbors on ideal lattice.
      typedef std::vector< std::vector< atat::rVector3d > > t_FirstNeighbors;
      t_FirstNeighbors fn;
      first_neighbors_( fn );
      for( size_t i(0); i < structure.lattice->sites.size(); ++i )
        std::sort( fn[i].begin(), fn[i].end(), details::Order() );

      centers.clear();
      
      // Creates a list of centers
      t_Atoms :: iterator i_atom = structure.atoms.begin();
      t_Atoms :: iterator i_atom_end = structure.atoms.end();

      for(types::t_unsigned index=0; i_atom != i_atom_end; ++i_atom, ++index )
        centers.push_back( AtomicCenter( structure, *i_atom, index ) );

      // builds an array of indices.
      std::vector< size_t > indices( centers.size() );
      {
        std::vector<size_t> :: iterator i_b = indices.begin();
        std::vector<size_t> :: iterator i_e = indices.end();
        for( size_t i(0); i_b != i_e; ++i_b, ++i ) *i_b = i;
      }
      t_Centers :: iterator i_center( centers.begin() );
      t_Centers :: iterator i_center_end = centers.end();
      for( i_center; i_center != i_center_end; ++i_center )
      {
        std::partial_sort( indices.begin(), 
                           indices.begin() + 5, 
                           indices.end(),
                           details::OrderI( i_center, centers ) );
        const site( structure.atoms[ indices.front() ].site );
        __ASSERT( site > 0 and site < structure.lattice->site.size(),
                  "sites not indexed.\n" )
        for( size_t i(1); i < 5; ++i )
        {
          t_Centers :: iterator i_bond( centers.begin() + indices[i] );
          i_center->bonds.push_back( t_Center::__make__iterator__( i_bond ) );
          const atat::rVector3d dfrac
          ( 
              inv_cell 
            * ( 
                  (const atat::rVector3d) *i_center 
                - (const atat::rVector3d) *i_bond
                + fn[site][i-1] 
              )
           ); 
          const atat::rVector3d frac
          (
            rint( dfrac(0) ),
            rint( dfrac(1) ),
            rint( dfrac(2) )
          );
          i_center->translations.push_back( frac );
          i_center->do_translates.push_back
          ( 
            atat::norm2(frac) > atat::zero_tolerance 
          );
        }
      }

      __DODEBUGCODE( check_tree(); )
      return true;
    } // Vff :: construct_bond_list

  } // namespace vff
} // namespace LaDa
