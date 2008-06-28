//
//  Version: $Id$
//

#include "gsl.h"

#include<boost/type_traits/is_same.hpp>

namespace Gsl
{
  Vector::Vector( std::vector<types::t_real> &_vec )
  {
    if( not boost::is_same< types::t_real, double>::value )
     __THROW_ERROR("Types are not equivalent, double != types::t_real.\n" )
    vector.size = _vec.size();
    vector.stride = 1;
    vector.data = &_vec[0];
    vector.block = &bloc;
    vector.owner = 0;
    bloc.size = vector.size;
    bloc.data = vector.data;
  }
  Matrix::Matrix( types::t_int _nrow, std::vector<types::t_real> &_mat )
  {
    if( not boost::is_same< types::t_real, double>::value )
     __THROW_ERROR("Types are not equivalent, double != types::t_real.\n" )
    bloc.size = _mat.size();
    bloc.data = &_mat[0];
    matrix.size1 = _nrow;
    matrix.size2 = _mat.size() / _nrow;
    __DOASSERT( matrix.size1 * matrix.size2 != _mat.size(),
                   "Matrix dimensions are not coherent: "
                << matrix.size1 << "x" << matrix.size2
                << "!=" << bloc.size << ".\n" )
    matrix.tda = matrix.size2;
    matrix.data = bloc.data;
    matrix.block = &bloc;
    matrix.owner = 0;
  }
}
