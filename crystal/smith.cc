//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#include <opt/smith_normal_form.h>
#include <opt/debug.h>

#include "smith.h"

namespace LaDa
{
  namespace Crystal
  { 
    atat::iVector3d get_smith_index( t_SmithTransform const &_transformation,
                                     atat::rVector3d  const &_pos )
    {
      namespace bt = boost::tuples;
      atat::iVector3d result;
      const atat::rVector3d pos( bt::get<0>( _transformation ) * _pos );
      const atat::iVector3d int_pos
      (
        types::t_int( rint( pos(0) ) ),
        types::t_int( rint( pos(1) ) ),
        types::t_int( rint( pos(2) ) )
      );
      for( size_t i(0); i < 3; ++i )
      {
        __DOASSERT
        (
          std::abs( pos(i) - types::t_real( int_pos(i) ) ) > 0.001, 
          "Structure is not ideal.\n"
        )
        result(i) = int_pos(i) % bt::get<1>(_transformation)(i);
        if( result(i) < 0 ) result(i) += bt::get<1>(_transformation)(i);
      }
      return result;
    }

    t_SmithTransform get_smith_transform( atat::rMatrix3d const &_lat_cell,
                                          atat::rMatrix3d const &_str_cell )
    {
      namespace bt = boost::tuples;
      t_SmithTransform result;
      atat::iMatrix3d left, right, smith;
      const atat::rMatrix3d inv_lat( !_lat_cell );
      const atat::rMatrix3d inv_lat_cell( inv_lat * _str_cell );
      atat::iMatrix3d int_cell;
      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
        {
          int_cell(i,j) = types::t_int( rint( inv_lat_cell(i,j) ) ); 
          __DOASSERT
          ( 
            std::abs( types::t_real( int_cell(i,j) ) - inv_lat_cell(i,j) ) > 0.01,
               "Input structure is not supercell of the lattice: \n" 
            << int_cell << "\n != \n" << inv_lat_cell << "\n\n"
          )
        }
      opt::smith_normal_form( smith, left, int_cell, right );
      for( size_t i(0); i < 3; ++i )
      {
        for( size_t j(0); j < 3; ++j )
          bt::get<0>( result )(i,j) = types::t_real( left(i,j) );
        bt::get<1>( result )(i) = smith(i,i);
      }
      bt::get<0>( result ) = bt::get<0>( result ) * ( !_lat_cell );

      return result;
    }

    void get_smith_map( Crystal::Structure const &_str, std::vector< std::vector<size_t> > &_out)
    {
      // finds map for atomic positions.
      Crystal::t_SmithTransform const transform( get_smith_transform(_str.lattice->cell, _str.cell) );

      LADA_ASSERT( _str.atoms.size() % _str.lattice->sites.size() == 0, "Unexpected number of atoms.\n" );
      size_t const N( _str.atoms.size() / _str.lattice->sites.size() );
      _out.clear();
      _out.resize(_str.lattice->sites.size(), std::vector<size_t>(N, N)); 

      Crystal::Structure::t_Atoms::const_iterator i_first = _str.atoms.begin();
      Crystal::Structure::t_Atoms::const_iterator const i_end = _str.atoms.end();
      for(size_t i(0); i_first != i_end; ++i_first, ++i)
      {
        LADA_ASSERT( i_first->site >= 0 and i_first->site < _str.lattice->sites.size(),
                     "Atoms in structure are not indexed with respect to sites.\n" );
        
        atat::rVector3d const pos( i_first->pos - _str.lattice->sites[i_first->site].pos );
        size_t const smith( Crystal::get_linear_smith_index(transform, pos) );
        LADA_ASSERT( smith < N, "Incoherent number of atoms per site." );
        LADA_ASSERT( _out[i_first->site][smith] == N, "Duplicate site." );
        _out[i_first->site][smith] = i;
      }
#     ifdef LADA_DEBUG
        for( size_t site(0); site < _str.lattice->sites.size(); ++site )
          for( size_t i(0); i < N; ++i )
            LADA_ASSERT(_out[site][i] != N, "Atomic site not mapped.")
#     endif
    }
  }
}
