//
//  Version: $Id: confsplit.cc 1207 2009-06-25 00:29:08Z davezac $
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<algorithm>
#include<boost/lambda/lambda.hpp>

#include<opt/fuzzy.h>

#include "neighbors.h"

namespace LaDa
{
  namespace Crystal
  {
    bool operator<( Neighbor const& _a, Neighbor const& _b )
      { return _a.distance < _b.distance; }

    void Neighbors :: create_neighbors_list_(t_Structure const& _structure)
    {
      const types::t_int N( _structure.atoms.size() );
      const types::t_int umax = nmax / _structure.atoms.size() + 1;
      neighbors_.clear();
      size_t size(0);
      
      atat::rMatrix3d const inv_cell( !_structure.cell );
      t_Structure::t_Atoms::const_iterator i_atom = _structure.atoms.begin();
      t_Structure::t_Atoms::const_iterator i_atom_end = _structure.atoms.end();
      Neighbor neighbor;
      neighbor.index = 0;
      for(; i_atom != i_atom_end; ++i_atom, ++neighbor.index ) 
      {
        atat::rVector3d const frac( inv_cell * (i_atom->pos - origin) );
        atat::rVector3d const centered
        ( 
          frac(0) - std::floor( frac(0) + 0.500000001e0 ),
          frac(1) - std::floor( frac(1) + 0.500000001e0 ),
          frac(2) - std::floor( frac(2) + 0.500000001e0 ) 
        );
        for( types::t_int x(-umax); x <= umax; ++x )
          for( types::t_int y(-umax); y <= umax; ++y )
            for( types::t_int z(-umax); z <= umax; ++z )
            {
               neighbor.pos = _structure.cell * ( frac + atat::rVector3d(x,y,z) );
               neighbor.distance = atat::norm( neighbor.pos );
               if( Fuzzy::is_zero( neighbor.distance ) ) continue;

               t_Neighbors :: iterator i_found 
               (
                 std::find_if
                 (
                   neighbors_.begin(), neighbors_.end(),
                   boost::lambda::constant(neighbor) < boost::lambda::_1  
                 ) 
               );

               
               if( i_found != neighbors_.end() ) 
               {
                 neighbors_.insert(i_found, neighbor);

                 if( size < nmax ) ++size;
                 else              neighbors_.pop_back();
                 continue;
               }
               
               if( size == nmax ) continue;

               ++size;
               neighbors_.push_back( neighbor );
            }
      } // loop over atoms.
    };

  } // namespace Crystal
} // namespace LaDa
