#include "LaDaConfig.h"

#include <iostream>
#include "gsl.h"

namespace LaDa
{
  namespace Gsl
  {
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
} // namespace LaDa
