//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>

#include "layerdepth.h"

namespace LaDa
{
  namespace Crystal 
  {

     void LayerDepth :: set( const Eigen::Matrix3d &_mat )
     {
       a0 = _mat.col(0).normalized();
       a1 = _mat.col(1).normalized();
       a2 = _mat.col(2).normalized();
       __DODEBUGCODE( isset = true; )
     } 

     void LayerDepth :: set( const Eigen::Vector3d& _vec )
     {
        a0 = vec.normalized();
        
        // First, lets create an orthonormal vector to _vec
        a1 = Fuzzy::eq( a0(0), 0e0 ) ? 
               ( Fuzzy::eq( a0(1), 0e0 ) ? 
                  Eigen::Vector3d(1, 0, 0):
                  Eigen::Vector3d(0, a0(2), -a0(1) )  ): 
               Eigen::Vector3d( -a0(2) -a0(1), a0(0), a0(0) );

        // Then, lets create another... 
        a2( 0 ) = a0(1) * a1(2) - a0(2) * a1(1);
        a2( 1 ) = a0(2) * a1(0) - a0(0) * a1(2);
        a2( 2 ) = a0(0) * a1(1) - a0(1) * a1(0);
     }
    
     //! Strict weak ordering operator.
     bool LayerDepth :: operator()( const Eigen::Vector3d& _first, 
                                    const Eigen::Vector3d& _second ) const
     {
       __ASSERT( not isset, "Crystal::LayerDepth has not been set.\n" )
       types::t_real a =  _first * a0; a -= std::floor( a - types::tolerance );
       types::t_real b =  _second * a0; b -= std::floor( b - types::tolerance );

       if ( not Fuzzy::eq( a, b ) ) return Fuzzy::le( a, b );
     
       a =   _first * a1; a -= std::floor( a - types::tolerance );
       b =   _second * a1; b -= std::floor( b - types::tolerance );
       if ( not Fuzzy::eq( a, b ) ) return Fuzzy::le( a, b );
       
       a =   _first * a2; a -= std::floor( a - types::tolerance );
       b =   _second * a2; b -= std::floor( b - types::tolerance );
       return Fuzzy::le( a, b );
     }
     //! Returns the depth.
     types::t_real LayerDepth :: operator()( const Eigen::Vector3d& _first ) const
     {
       __ASSERT( not isset, "Crystal::LayerDepth has not been set.\n" )
       const types::t_real a =  _first * a0; 
       return a - std::floor( a - types::tolerance );
     }

  } // namespace Crystal
} // namespace LaDa
