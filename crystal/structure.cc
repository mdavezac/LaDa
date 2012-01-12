#include "LaDaConfig.h"

#include <algorithm>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple_io.hpp>

#include <Eigen/LU>

#include <opt/tinyxml.h>
#include <opt/ndim_iterator.h>
#include <physics/physics.h>
#include <math/misc.h>

#include "structure.h"
#include "fill_structure.h"
#include "epi_structure.h"
#include "fourier.h"
#include "smith.h"

namespace LaDa
{

  namespace crystal 
  {
    // Transforms a structure according to an affine transformation.
    Structure Structure::transform(math::Affine3d const &_affine) const
    {
      Structure result = copy();
      result.cell() = _affine.linear() * cell();
      iterator i_first = result.begin();
      iterator const i_end = result.end();
      for(; i_first != i_end;  ++i_first)
        i_first->pos() = _affine * i_first->pos();
      return result;
    }
  } // namespace Crystal

} // namespace LaDa
