#ifndef LADA_MATH_COMPARE_NORMS_H
#define LADA_MATH_COMPARE_NORMS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "fuzzy.h"

namespace LaDa
{
  namespace math
  { 
    //! Functor for sorting atat real vectors
    class CompareNorms
    {
      //! Vector to which to compare
      const math::rVector3d vec ;
      public:
        //! Constructor and initializer
        norm_compare( math::rVector3d &_vec ) : vec(_vec) {};
        //! Constructor and initializer
        norm_compare() : vec(0,0,0) {};
        //! Copy Constructor
        norm_compare( const norm_compare &_c ) : vec(_c.vec) {};
        //! \brief returns true if the norm of \a _a - norm::compare::vec is
        //!        smaller than the norm pf \a _b - norm::compare::vec.
        bool operator()( math::rVector3d &_a, math::rVector3d &_b ) const
          { return math::le( (_a - vec).squaredNorm(), (_b - vec).squaredNorm() ); }
        //! Returns true if \a _a = norm::compare::vec.
        bool operator()( math::rVector3d &_a ) const
          { return math::is_zero( (_a - vec).squaredNorm() ); } 
    };
  }
} // namespace LaDa
#endif
