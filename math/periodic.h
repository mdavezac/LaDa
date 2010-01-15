#ifndef LADA_MATH_PERIODIC_H
#define LADA_MATH_PERIODIC_H

#include <Eigen/Core>

#include "is_integer.h"

namespace LaDa
{
  namespace math
  {
    inline bool are_periodic_images( Eigen::Vector3d const &_a,
                                     Eigen::Vector3d const &_b, 
                                     Eigen::Matrix3d const &_inv_cell )
    {
      return is_integer( inv_cell * (_a - _b) );
    }
  }
}

#endif
