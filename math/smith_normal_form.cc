#include "LaDaConfig.h"

#include <opt/types.h>
#include "smith_normal_form.h"
#include "fuzzy.h"
#include "exceptions.h"

namespace LaDa
{

  namespace math
  {

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
      result.get<0>() = left.cast<types::t_real>() * (!_unitcell);
      result.get<1>() = smith.diagonal();

      return result;
    }

    template<class T>
      bool one_nonzero(Eigen::MatrixBase<T> const &_in)
      {
        for(size_t i(0), n(0); i < _in.size(); ++i)
          if(_in(i) != 0) 
          {
            if(n == 1) return false;
            n = 1;
          }
        return true;
      };
    template<class T>
      bool check_nonzero(Eigen::MatrixBase<T> const &_in, size_t _index)
      {
        size_t n(0);
        for(size_t i(0); i < _in.rows(); ++i)
          if(_in(i, _index) != 0) 
          {
            if(n == 1) return false;
            n = 1;
          }
        for(size_t i(0); i < _in.cols(); ++i)
          if(i != _index and _in(_index, i) != 0) 
          {
            if(n == 1) return false;
            n = 1;
          }
        return true;
      };

    template<class T>
      void get_min_max(Eigen::MatrixBase<T> const &_in, types::t_int &_max, types::t_int &_min)
      {
        _max = 0;
        for(size_t k(1); k < _in.size(); ++k)
          if(std::abs(_in(k)) > std::abs(_in(_max))) _max = k;
        types::t_int k(_in.size()-1);
        for(; k >= 0 and _in(k) == 0; --k);
        _min = k;
        for(--k; k >= 0; --k)
          if(_in(k) != 0 and std::abs(_in(k)) < std::abs(_in(_min))) _min = k;

      }

    template<class T>
      void check_same_magnitude( Eigen::MatrixBase<T> const &_in,
                                 types::t_int &_max, types::t_int &_min, size_t _index )
      {
        if(std::abs(_in(_min, _index)) != std::abs(_in(_max, _index))) return;
        size_t n0(0), n1(0);
        for(size_t i(0); i < _in.rows(); ++i)
        {
          if(_in(_max, i)) ++n0;
          if(_in(_min, i)) ++n1;
        }
        if(n0 < n1 or (n0 == n1 and _in(_max, _index) < _in(_max, _index)) ) std::swap(_min, _max);
      }

    //! Column 0 to have only one non-zero positive component placed at origin.
    template<class T0, class T1, class T2, class T3>
      void smith_col_impl( Eigen::MatrixBase<T0> const &_out, 
                           Eigen::MatrixBase<T1> &_left,
                           Eigen::MatrixBase<T2> &_smith, 
                           Eigen::MatrixBase<T3> &_right, 
                           size_t _index )
      {
        while(not one_nonzero(_smith.col(_index)))
        {
          // find min/max elements.
          types::t_int maxelem, minelem;
          get_min_max(_smith.col(_index), maxelem, minelem); 
          check_same_magnitude(_smith, maxelem, minelem, _index);
          // Remove multiple from column.
          types::t_int const multiple = _smith(maxelem, _index) / _smith(minelem, _index);
          _smith.row(maxelem) -= multiple * _smith.row(minelem);
          _left.row(maxelem) -= multiple * _left.row(minelem);
        }
        if(_smith(_index, _index) == 0) 
        {
          size_t k(0);
          for(; k < _smith.rows() and _smith(k, _index) == 0; ++k);
          if(k == _smith.rows()) BOOST_THROW_EXCEPTION(error::internal());
          _smith.row(k).swap(_smith.row(_index));
          _left.row(k).swap(_left.row(_index));
        }
        if(_smith(_index, _index) < 0)
        {
          _smith.row(_index) *= -1;
          _left.row(_index) *= -1;
        }
      }

    //! Row 0 to have one only non-zero and positive component placed at origin.
    template<class T0, class T1, class T2, class T3>
      void smith_row_impl( Eigen::MatrixBase<T0> const &_out, 
                           Eigen::MatrixBase<T1> &_left,
                           Eigen::MatrixBase<T2> &_smith,
                           Eigen::MatrixBase<T3> &_right,
                           size_t _index )
      {
        while(not one_nonzero(_smith.row(_index)))
        {
          // find min/max elements.
          types::t_int maxelem, minelem;
          get_min_max(_smith.row(_index).transpose(), maxelem, minelem); 
          check_same_magnitude(_smith.transpose(), maxelem, minelem, _index);
          // Remove multiple from column.
          types::t_int const multiple = _smith(_index, maxelem) / _smith(_index, minelem);
          _smith.col(maxelem) -= multiple * _smith.col(minelem);
          _right.col(maxelem) -= multiple * _right.col(minelem);
        }
        if(_smith(_index, _index) == 0) 
        {
          size_t k(0);
          for(; k < _smith.cols() and _smith(_index, k) == 0; ++k);
          if(k == _smith.cols()) BOOST_THROW_EXCEPTION(error::internal());
          _smith.col(k).swap(_smith.col(_index));
          _right.col(k).swap(_right.col(_index));
        }
        if(_smith(_index, _index) < 0)
        {
          _smith.col(_index) *= -1;
          _right.col(_index) *= -1;
        }
      }


    //! Makes matrix diagonal. Does not order diagonal values correctly yet.
    template<class T0, class T1, class T2, class T3>
      void smith_impl_( Eigen::MatrixBase<T0> const &_out, 
                        Eigen::MatrixBase<T1> &_left,
                        Eigen::MatrixBase<T2> &_right,
                        Eigen::MatrixBase<T3> &_smith )
      {
        size_t const nrows = _smith.rows();
        size_t const ncols = _smith.cols();
        Eigen::Matrix<typename Eigen::MatrixBase<T0>::Scalar, Eigen::Dynamic, Eigen::Dynamic> old(nrows, ncols);
        for(size_t index(0); index < _smith.rows()-1; ++index)
          do
          {
            smith_col_impl(_out, _left, _smith, _right, index);
            smith_row_impl(_out, _left, _smith, _right, index);

            size_t maxrow = _smith.rows();
            types::t_int const diag = _smith(index, index);
            types::t_int maxmod = 0;
            for(size_t i(index+1); i < _smith.rows(); ++i)
              for(size_t j(index+1); j < _smith.cols(); ++j)
              {
                if(_smith(i, j) % diag == 0) continue;
                else if(maxmod == 0) { maxrow = i; maxmod = std::abs(_smith(i,j) % diag); }
                else if(std::abs(_smith(i,j) % diag) > maxmod)
                  { maxrow = i; maxmod = std::abs(_smith(i,j) % diag); }
              }
            if(maxmod != 0) 
            {
              _smith.row(index) += _smith.row(maxrow);
              _left.row(index)  += _left.row(maxrow);
            }
          } while( not check_nonzero(_smith, index) );
        if(_smith(_smith.rows()-1, _smith.cols()-1) < 0)
        {
          _smith.row(_smith.rows()-1) *= -1;
          _left.row(_smith.rows()-1) *= -1;
        }
      }
    
   void smith_normal_form( iMatrix3d& _S, iMatrix3d & _L,
                           const iMatrix3d& _M, iMatrix3d &_R )
   {
     // set up the system _out = _left * _smith * _right.
     _L  = iMatrix3d::Identity();
     _R = iMatrix3d::Identity();
     _S = _M;
     smith_impl_(_M, _L, _R, _S);
   }
  } // namespace Crystal

} // namespace LaDa
