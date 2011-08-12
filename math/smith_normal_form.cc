#include "LaDaConfig.h"
#include "FCMangle.h"

#include <opt/types.h>
#include "smith_normal_form.h"
#include "fuzzy.h"
#include "exceptions.h"

//! \cond
extern "C" void FC_GLOBAL(smithnormalform, SMITHNORMALFORM)
                         ( const int*, int*, int*, int* );
//! \endcond

namespace LaDa
{

  namespace math
  {
    
    void smith_normal_form( math::iMatrix3d& _S, math::iMatrix3d & _L,
                            const math::iMatrix3d& _M, math::iMatrix3d &_R )
    {
      types::t_int s[9], l[9], m[9], r[9];
      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
        {
          m[ j * 3 + i ] = _M(i,j);
          s[ j * 3 + i ] = 0;
          l[ j * 3 + i ] = 0;
          r[ j * 3 + i ] = 0;
        }
      FC_GLOBAL(smithnormalform, SMITHNORMALFORM)( m, l, s, r );
      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
        {
          _S(i,j) = s[ j * 3 + i ];
          _L(i,j) = l[ j * 3 + i ];
          _R(i,j) = r[ j * 3 + i ];
        }
    }

    math::iVector3d smith_index( t_SmithTransform const &_transformation,
                                 math::rVector3d  const &_pos )
    {
      namespace bt = boost::tuples;
      math::iVector3d result;
      const math::rVector3d pos( bt::get<0>( _transformation ) * _pos );
      const math::iVector3d int_pos
      (
        types::t_int( rint( pos(0) ) ),
        types::t_int( rint( pos(1) ) ),
        types::t_int( rint( pos(2) ) )
      );
      for( size_t i(0); i < 3; ++i )
      {
#       ifdef LADA_DEBUG
          if( math::neq(pos(i), types::t_real(int_pos(i))) )
            BOOST_THROW_EXCEPTION(error::off_lattice_position());
#       endif
        result(i) = int_pos(i) % bt::get<1>(_transformation)(i);
        if( result(i) < 0 ) result(i) += bt::get<1>(_transformation)(i);
      }
      return result;
    }

    t_SmithTransform smith_transform( math::rMatrix3d const &_unitcell,
                                      math::rMatrix3d const &_supercell )
    {
      namespace bt = boost::tuples;
      t_SmithTransform result;
      math::iMatrix3d left, right, smith;
      const math::rMatrix3d inv_lat( !_unitcell );
      const math::rMatrix3d inv_lat_cell( inv_lat * _supercell );
      math::iMatrix3d int_cell;
      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
        {
          int_cell(i,j) = types::t_int( rint( inv_lat_cell(i,j) ) ); 
          if( math::neq(types::t_real(int_cell(i,j)), inv_lat_cell(i,j), 1e-2) )
            BOOST_THROW_EXCEPTION( error::not_a_supercell() );
        }
      math::smith_normal_form( smith, left, int_cell, right );
      for( size_t i(0); i < 3; ++i )
      {
        for( size_t j(0); j < 3; ++j )
          bt::get<0>( result )(i,j) = types::t_real( left(i,j) );
        bt::get<1>( result )(i) = smith(i,i);
      }
      bt::get<0>( result ) = bt::get<0>( result ) * ( !_unitcell );

      return result;
    }
   
  } // namespace Crystal

} // namespace LaDa
