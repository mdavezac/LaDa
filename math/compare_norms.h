#ifndef LADA_MATH_COMPARE_NORMS_H
#define LADA_MATH_COMPARE_NORMS_H

#include "LaDaConfig.h"

#include "fuzzy.h"

namespace LaDa
{
  namespace math
  { 
    //! Functor for sorting eigen real vectors
    class CompareNorms
    {
      //! Vector to which to compare
      const math::rVector3d vec ;
      public:
        //! Constructor and initializer
        CompareNorms( math::rVector3d const &_vec ) : vec(_vec) {};
        //! Constructor and initializer
        CompareNorms() : vec(0,0,0) {};
        //! Copy Constructor
        CompareNorms( const CompareNorms &_c ) : vec(_c.vec) {};
        //! \brief returns true if the norm of \a _a - norm::compare::vec is
        //!        smaller than the norm pf \a _b - norm::compare::vec.
        bool operator()( math::rVector3d &_a, math::rVector3d &_b ) const
          { return math::le( (_a - vec).squaredNorm(), (_b - vec).squaredNorm() ); }
        //! Returns true if \a _a = norm::compare::vec.
        bool operator()( math::rVector3d &_a ) const
          { return math::is_null( (_a - vec).squaredNorm() ); } 
    };
  }
} // namespace LaDa
#endif
