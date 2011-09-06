#include "LaDaConfig.h"

#include <limits>
#include <set>

#include <math/misc.h>

#include "periodic_dnc.h"

namespace LaDa 
{
  namespace crystal 
  {
      math::iVector3d DnCBoxes::guess_mesh( math::rMatrix3d const &_cell, size_t _N, size_t _nperbox ) const
      {
        types::t_real const alpha = 0.5;
        types::t_real const c0 = std::sqrt( _cell.col(0).squaredNorm() );
        types::t_real const c1 = std::sqrt( _cell.col(1).squaredNorm() );
        types::t_real const c2 = std::sqrt( _cell.col(2).squaredNorm() );

        int const Nboxes = int(std::floor(types::t_real(_N)/types::t_real(_nperbox) + 0.5));
        if(Nboxes == 0) return math::iVector3d::Ones();

        math::iVector3d result(1,1,1);
        types::t_real mini = _N / _nperbox * (c0+c1+c2);
        for(size_t n0(1); n0 <= Nboxes; ++n0)
          for(size_t n1(1); n1 <= Nboxes; ++n1)
          {
            size_t n2 = Nboxes / n0 / n1;
            if(n2 == 0) continue;

            types::t_real const a =
                 std::abs(c0/types::t_real(n0) - c1/types::t_real(n1))
               + std::abs(c0/types::t_real(n0) - c2/types::t_real(n2))
               + std::abs(c1/types::t_real(n1) - c2/types::t_real(n2));
            if(a < mini) 
            {
              result(0) = n0;
              result(1) = n1;
              result(2) = n2;
              mini = a;
            }
          }
        return result;
      }

  } // namespace Crystal

} // namespace LaDa

