#ifndef LADA_MATH_EIGENOPS_H
#define LADA_MATH_EIGENOPS_H

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Geometry>

#include <python/types.h>

namespace LaDa
{
  namespace math
  {
    //! \f$\pi\f$
    const types::t_real pi = 3.1415926535897932384626433832795028841971693993751058209749445920;

    //! 3d-vector of reals. 
    typedef Eigen::Matrix<types::t_real, 3, 1> rVector3d;
    //! 3d-vector of integers. 
    typedef Eigen::Matrix<types::t_int, 3, 1> iVector3d;
    //! 6d-vector of reals. 
    typedef Eigen::Matrix<types::t_real, 4, 1> rVector4d;
    //! 6d-vector of reals. 
    typedef Eigen::Matrix<types::t_real, 5, 1> rVector5d;
    //! 6d-vector of reals. 
    typedef Eigen::Matrix<types::t_real, 6, 1> rVector6d;

    //! 3d-matrix of reals. 
    typedef Eigen::Matrix<types::t_real, 3, 3> rMatrix3d;
    //! 3d-vector of integers. 
    typedef Eigen::Matrix<types::t_int, 3, 3> iMatrix3d;

    //! \typedef type of the angle axis object to initialize roations.
    typedef Eigen::AngleAxis<types::t_real> AngleAxis;
    //! \typedef type of the translation objects.
    typedef Eigen::Translation<types::t_real, 3> Translation;
#   ifndef LADA_WITH_EIGEN3 
      //! \typedef type of the affine transformations.
      typedef Eigen::Transform<types::t_real, 3> Affine3d;
#   else
      //! \typedef type of the affine transformations.
      typedef Eigen::Transform<types::t_real, 3, Eigen::Isometry> Affine3d;
#   endif
  } // namespace math
} // namespace LaDa


namespace Eigen
{
  //! Real type.
  typedef LaDa::types::t_real t_real;
  //! Integer type.
  typedef LaDa::types::t_int t_int;

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

# ifndef LADA_WITH_EIGEN3 
    //! Transpose operation of real matrix.
    inline Transpose< Matrix<t_real, 3, 3> > operator~(Matrix<t_real, 3, 3> const &_mat)
      { return _mat.transpose(); }
    //! Transpose operation of integer matrix.
    inline Transpose< Matrix<t_int, 3, 3> > operator~(Matrix<t_int, 3, 3> const &_mat)
      { return _mat.transpose(); }
# else
    //! Transpose operation of matrix.
    template<class T_DERIVED>
      inline typename MatrixBase<T_DERIVED>::ConstTransposeReturnType
        operator~(MatrixBase<T_DERIVED> const &_mat) { return _mat.transpose(); }
# endif
}
#endif
