#include "LaDaConfig.h"
#include "FCMangle.h"

#include <boost/math/common_factor_rt.hpp>
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
    
    void smith_normal_form( iMatrix3d& _S, iMatrix3d & _L,
                            const iMatrix3d& _M, iMatrix3d &_R )
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

    iVector3d smith_index( t_SmithTransform const &_transformation,
                           rVector3d  const &_pos )
    {
      namespace bt = boost::tuples;
      iVector3d result;
      const rVector3d pos( bt::get<0>( _transformation ) * _pos );
      const iVector3d int_pos
      (
        types::t_int( rint( pos(0) ) ),
        types::t_int( rint( pos(1) ) ),
        types::t_int( rint( pos(2) ) )
      );
      for( size_t i(0); i < 3; ++i )
      {
#       ifdef LADA_DEBUG
          if( neq(pos(i), types::t_real(int_pos(i))) )
            BOOST_THROW_EXCEPTION(error::off_lattice_position());
#       endif
        result(i) = int_pos(i) % bt::get<1>(_transformation)(i);
        if( result(i) < 0 ) result(i) += bt::get<1>(_transformation)(i);
      }
      return result;
    }

    t_SmithTransform smith_transform( rMatrix3d const &_unitcell,
                                      rMatrix3d const &_supercell )
    {
      if(std::abs(_unitcell.determinant()) < 1e-8)
        BOOST_THROW_EXCEPTION(error::singular_matrix() << error::string("Unit-cell is singular."));
      if(std::abs(_supercell.determinant()) < 1e-8)
        BOOST_THROW_EXCEPTION(error::singular_matrix() << error::string("Supercell is singular."));
      if(_unitcell.determinant() < 1e-8)
        BOOST_THROW_EXCEPTION( error::negative_volume()
                                 << error::string("Unit-cell volume is negative."));
      namespace bt = boost::tuples;
      if(_supercell.determinant() < 1e-8)
        BOOST_THROW_EXCEPTION( error::negative_volume()
                                 << error::string("Supercell volume is negative."));
      namespace bt = boost::tuples;
      t_SmithTransform result;
      iMatrix3d left, right, smith;
      const rMatrix3d inv_lat( !_unitcell );
      const rMatrix3d inv_lat_cell( inv_lat * _supercell );
      iMatrix3d int_cell;
      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
        {
          int_cell(i,j) = types::t_int( rint( inv_lat_cell(i,j) ) ); 
          if( neq(types::t_real(int_cell(i,j)), inv_lat_cell(i,j), 1e-2) )
            BOOST_THROW_EXCEPTION( error::not_a_supercell() );
        }
      smith_normal_form( smith, left, int_cell, right );
      for( size_t i(0); i < 3; ++i )
      {
        for( size_t j(0); j < 3; ++j )
          bt::get<0>( result )(i,j) = types::t_real( left(i,j) );
        bt::get<1>( result )(i) = smith(i,i);
      }
      bt::get<0>( result ) = bt::get<0>( result ) * ( !_unitcell );

      return result;
    }

    template<class T0, class T1, class T2, class T3>
    void smith_transform( Eigen::MatrixBase<T0> const &_out, 
                          Eigen::MatrixBase<T1> &_left,
                          Eigen::MatrixBase<T2> &_right,
                          Eigen::MatrixBase<T2> &_smith )
    {
      // loop over diagonal elements.
      size_t jt = 0;
      { // step 1 in wikipedia page, for columns.
        
        // find column index with at least one non-zero entry.
        for(; jt < M.cols() and M.col(jt) != 0; ++jt);
        if(jt == M.cols()) return;
        // swap two rows if condition is true.
        if(M(t, jt) == 0) 
        {
          size_t k(1);
          for(; k < M.rows() and M(k, jt) == 0; ++k);
#         ifdef LADA_DEBUG
            if(k == M.rows())
              BOOST_THROW_EXCEPTION( error::internal()
                                       << error::string("Could not find non-zero element."));
#         endif
          M.row(0).swap(M.row(k));
          _left.row(0).swap(_left.row(k));
        }
      }

      { // step 2 in wikipedia page, for columns.
        size_t k(1);
        for(; k < M.rows() and M(t, jt) % M(k, jt) == 0; ++k);
        if(k != M.rows())
        {
        }
      }

    }

    void smith_transform( iMatrix3d const &_out, iMatrix3d &_left,
                          iMatrix3d const &_right, iMatrix3d &_smith )
    {
      // set up the system _out = _left * _smith * _right.
      _left  = iMatrix3d::Identity();
      _right = iMatrix3d::Identity();
      _smith = _out;

      smith_transform_impl(_out, _left, _right, _smith);
    }

   
  } // namespace Crystal

} // namespace LaDa
