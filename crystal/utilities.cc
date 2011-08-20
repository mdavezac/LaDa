#include "LaDaConfig.h"

#include <math/misc.h>
#include <limits>

#include <math/fuzzy.h>
#include "utilities.h"


namespace LaDa
{
  namespace crystal 
  {
    math::rVector3d into_cell( math::rVector3d const &_vec, 
                               math::rMatrix3d const &_cell, 
                               math::rMatrix3d const &_inv)
    {
      math::rVector3d result( _inv * _vec );
      result(0) -= std::floor(result(0)+types::roundoff);
      result(1) -= std::floor(result(1)+types::roundoff);
      result(2) -= std::floor(result(2)+types::roundoff);
      return _cell * result;
    }

    math::rVector3d zero_centered( math::rVector3d const &_vec, 
                                   math::rMatrix3d const &_cell, 
                                   math::rMatrix3d const &_inv)
    {
      math::rVector3d result( _inv * _vec );
      result(0) -= std::floor(5e-1+result(0)+types::roundoff);
      result(1) -= std::floor(5e-1+result(1)+types::roundoff);
      result(2) -= std::floor(5e-1+result(2)+types::roundoff);
      // numerical stability check.
      if( math::eq(result(0), 5e-1) ) result(0) = -5e-1;
      else if( math::lt(result(0), -5e-1)) result(0) += 1e0;
      if( math::eq(result(1), 5e-1) ) result(1) = -5e-1;
      else if( math::lt(result(1), -5e-1)) result(1) += 1e0;
      if( math::eq(result(2), 5e-1) ) result(2) = -5e-1;
      else if( math::lt(result(2), -5e-1)) result(2) += 1e0;
      return _cell * result;
    }

    math::rVector3d into_voronoi( math::rVector3d const &_vec, 
                                  math::rMatrix3d const &_cell, 
                                  math::rMatrix3d const &_inv)
    {
      math::rVector3d result( _inv * _vec );
      result(0) -= std::floor(result(0)+types::roundoff);
      result(1) -= std::floor(result(1)+types::roundoff);
      result(2) -= std::floor(result(2)+types::roundoff);
      // numerical stability check.
      if( math::eq(result(0), -1e0) or math::eq(result(0), 1e0) ) result(0) = 0e0;
      if( math::eq(result(1), -1e0) or math::eq(result(1), 1e0) ) result(1) = 0e0;
      if( math::eq(result(2), -1e0) or math::eq(result(2), 1e0) ) result(2) = 0e0;
      math::rVector3d const orig(result);
      types::t_real min_norm = (_cell*orig).squaredNorm();
      for(int i(-1); i < 2; ++i)
        for(int j(-1); j < 2; ++j)
          for(int k(-1); k < 2; ++k)
          {
            math::rVector3d const translated = orig + math::rVector3d(i,j,k);
            types::t_real const d( (_cell*translated).squaredNorm() );
            if( math::gt(min_norm, d) )
            {
              min_norm = d;
              result = translated;
            }
          }
      return _cell * result;
    }

  } // namespace Crystal
} // namespace LaDa
