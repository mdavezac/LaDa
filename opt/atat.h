//
//  Version: $Id$
//
#ifndef __ATAT_H_
#define __ATAT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <atat/vectmac.h>

#include "fuzzy.h"

namespace LaDa
{
  namespace atat
  { 
    //! Functor for sorting atat real vectors
    class norm_compare
    {
      //! Vector to which to compare
      const rVector3d vec ;
      public:
        //! Constructor and initializer
        norm_compare( const rVector3d &_vec ) : vec(_vec) {};
        //! Constructor and initializer
        norm_compare() : vec(0,0,0) {};
        //! Copy Constructor
        norm_compare( const norm_compare &_c ) : vec(_c.vec) {};
        //! \brief returns true if the norm of \a _a - norm::compare::vec is
        //!        smaller than the norm pf \a _b - norm::compare::vec.
        bool operator()( const rVector3d &_a, const rVector3d &_b ) const
          { return Fuzzy::le( norm2( _a - vec ), norm2( _b - vec ) ); }
        //! Returns true if \a _a = norm::compare::vec.
        bool operator()( const rVector3d &_a ) const
          { return Fuzzy::eq( norm2( _a - vec ), 0.0 ); } 
    };
  }
} // namespace LaDa
#endif
