//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>
#include <math/eigen.h>

#include "layerdepth.h"

namespace LaDa
{
  namespace Crystal 
  {
     //! Strict weak ordering operator.
     bool LayerDepth :: operator()( const math::rVector3d& _first, 
                                    const math::rVector3d& _second ) const
     {
       math::rVector3d a = matrix_.inverse() * _first;
       math::rVector3d b = matrix_.inverse() * _second;
       a(0) = a(0) - std::floor(a(0) + types::tolerance);
       if( math::eq( a(0), 1e0 ) ) a(0) -= 1e0;
       b(0) = b(0) - std::floor(b(0) + types::tolerance);
       if( math::eq( b(0), 1e0 ) ) b(0) -= 1e0;

       if ( not math::eq( a(0), b(0) ) ) return math::gt( a(0), b(0) );
       if ( not math::eq( a(1), b(1) ) ) return math::gt( a(1), b(1) );
       if ( not math::eq( a(2), b(2) ) ) return math::gt( a(2), b(2) );
       return false;
     }
     //! Returns the depth.
     types::t_real LayerDepth :: operator()( const math::rVector3d& _first ) const
     {
       math::rVector3d a = matrix_.inverse() * _first;
       a(0) = a(0) - std::floor(a(0) + types::tolerance);
       return ( math::eq( a(0), 1e0 ) ) ? a(0) - 1e0: a(0);
     }

  } // namespace Crystal
} // namespace LaDa
