#ifndef LADA_MATH_EIGENOPS_H
#define LADA_MATH_EIGENOPS_H

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Geometry>

#include <opt/types.h>

namespace LaDa
{
  namespace math
  {
    //! 3d-vector of reals. 
    typedef Eigen::Matrix<types::t_real, 3, 1> rVector3d;
    //! 3d-vector of integers. 
    typedef Eigen::Matrix<types::t_int, 3, 1> iVector3d;

    //! 3d-matrix of reals. 
    typedef Eigen::Matrix<types::t_real, 3, 3> rMatrix3d;
    //! 3d-vector of integers. 
    typedef Eigen::Matrix<types::t_int, 3, 3> iMatrix3d;
  } // namespace math
} // namespace LaDa


namespace Eigen
{
  //! Real type.
  typedef LaDa::types::t_real t_real;
  //! Integer type.
  typedef LaDa::types::t_int t_int;

// //! Dot product of real vectors.
// inline t_real operator*(Matrix<t_real, 3, 1> const &_a, Matrix<t_real, 3, 1> const &_b)
//   { return _a.dot(_b); }
// //! Dot product of integer vectors.
// inline t_real operator*(Matrix<t_int, 3, 1> const &_a, Matrix<t_int, 3, 1> const &_b)
//   { return _a.dot(_b); }

  //! Cross product of real vectors.
  inline Matrix<t_real, 3, 1> operator^(Matrix<t_real, 3, 1> const &_a, Matrix<t_real, 3, 1> const &_b)
    { return _a.cross(_b); }
  //! Cross product of integer vectors.
  inline Matrix<t_int, 3, 1> operator^(Matrix<t_int, 3, 1> const &_a, Matrix<t_int, 3, 1> const &_b)
    { return _a.cross(_b); }

  //! \brief Inverse operation of real matrix.
  //! \note Probably slower than using eigen because of return type.
  inline Matrix<t_real, 3, 3> operator!(Matrix<t_real, 3, 3> const &_mat)
    { return _mat.inverse(); }

  //! Transpose operation of real matrix.
  inline Eigen::Transpose< Matrix<t_real, 3, 3> > operator~(Matrix<t_real, 3, 3> const &_mat)
    { return _mat.transpose(); }
  //! Transpose operation of integer matrix.
  inline Eigen::Transpose< Matrix<t_int, 3, 3> > operator~(Matrix<t_int, 3, 3> const &_mat)
    { return _mat.transpose(); }
}
#endif
