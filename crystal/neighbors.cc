#include "LaDaConfig.h"

#include <opt/debug.h>
#include <math/fuzzy.h>

#include "neighbors.h"

namespace LaDa
{

  namespace Crystal 
  {

    //! \brief Strict weak ordering for distance from origin.
    //! \details Comparison with origin is always false.
    //!          This way, ends up a the end when sorting.
    struct Order
    {
      Order( types::t_real r ) : m(r) {}
      Order( const Order &_c ) : m(_c.m) {}
      bool operator()( const math::rVector3d& _a, const math::rVector3d &_b )
      {
        const types::t_real a = _a.squaredNorm();
        const types::t_real b = _b.squaredNorm();
        if( std::abs( a - b ) > m ) return a < b;
        if( math::neq( _a[0], _b[0] ) ) return _a[0] < _b[0];
        if( math::neq( _a[1], _b[1] ) ) return _a[1] < _b[1];
        return _a[2] < _b[2];
      }
      types::t_real m;
    };

    void find_first_neighbors( std::vector< math::rVector3d > &_positions,
                               const math::rMatrix3d &_cell,
                               const size_t _n )
    {
      LADA_NASSERT( _positions.size() == 0, "Position vector is empty.\n" )
      // first recenters around first atom.
      const math::rMatrix3d inv_cell( _cell.inverse() );
      const math::rVector3d origin( _positions[0] );
      types::t_real mindist = -1e0;
      foreach( math::rVector3d &pos, _positions )
      {
        pos -= origin;
        pos = inv_cell * pos;
        for( size_t i(0); i < 3; ++i )
          pos[i] = pos[i] - rint(pos[i]);
        pos = _cell * pos;
        const types::t_real d(pos.squaredNorm());
        if( d < types::tolerance ) continue;
        if( mindist > d  or mindist < 0e0 ) mindist = d;
      }
      math::iVector3d range
      (
        std::sqrt( _n * mindist / _cell.col(0).squaredNorm() ),
        std::sqrt( _n * mindist / _cell.col(1).squaredNorm() ),
        std::sqrt( _n * mindist / _cell.col(2).squaredNorm() )
      ); 
      _positions.reserve( _n + (2*range(0)+2) * (2*range(1)+1) * (2*range(2)+1) );
      {
        const std::vector< math::rVector3d > copy( _positions );
        std::vector< math::rVector3d > :: const_iterator i_pos = copy.begin();
        std::vector< math::rVector3d > :: const_iterator i_pos_end = copy.end();
        for(; i_pos != i_pos_end; ++i_pos )
        {
          for( types::t_int i(-range(0) ); i <= range(0); ++i )
            for( types::t_int j(-range(1) ); j <= range(1); ++j )
              for( types::t_int k(-range(2) ); k <= range(2); ++k )
              {
                if( i == 0 and j == 0 and k == 0 ) continue;
                const math::rVector3d d( (*i_pos) + _cell * math::rVector3d( i,j,k )  );
                _positions.push_back( d );
              }
        }
      }

      // then resizes list of atoms to _n 
      std::partial_sort( _positions.begin(), _positions.begin() + _n + 1, 
                         _positions.end(), Order( mindist * 0.25 ) );
      std::swap( _positions.front(), _positions[_n] );
      _positions.resize( _n );
    }

  } // namespace Crystal

} // namespace LaDa
