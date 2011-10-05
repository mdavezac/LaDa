#ifndef LADA_MATH_GRUBER_H
#define LADA_MATH_GRUBER_H

#include "LaDaConfig.h"

#include "eigen.h"

namespace LaDa
{
  namespace math
  {
    //! \brief Computes the Gruber cell of a lattice.
    //! \details Shamelessly adapted from the computational crystallography
    //!          tool box.
    rMatrix3d gruber( rMatrix3d const &_in, size_t itermax = 0,
                      types::t_real _tol = types::tolerance );
    //! \brief Computes the parameters of the metrical matrix associated with
    //!        the gruber cell.
    Eigen::Matrix<types::t_real, 6, 1>
      gruber_parameters( rMatrix3d const &_in, size_t itermax = 0,
                         types::t_real _tol = types::tolerance);
    namespace details
    {
      //! \brief function which needs to be compiled without optimization.
      //! \throw error::singular_matrix if the matrix is singular.
      bool no_opt_change_test(types::t_real _new, types::t_real _last);
    }
  }
}

#endif

