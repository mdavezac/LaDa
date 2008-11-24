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

     void LayerDepth :: set( const atat::rMatrix3d &_mat )
     {
       a0 = 1e0 / atat::norm2( _mat.get_column(0) ) * _mat.get_column(0);
       a1 = 1e0 / atat::norm2( _mat.get_column(1) ) * _mat.get_column(1);
       a2 = 1e0 / atat::norm2( _mat.get_column(2) ) * _mat.get_column(2);
       __DODEBUGCODE( isset = true; )
     } 

     void LayerDepth :: set( const atat::rVector3d& _vec )
     {
        a0 = 1e0 / atat::norm2( _vec ) * _vec ;
        
        // First, lets create an orthonormal vector to _vec
        a1 = Fuzzy::eq( a0(0), 0e0 ) ? 
               ( Fuzzy::eq( a0(1), 0e0 ) ? 
                  atat::rVector3d(1, 0, 0):
                  atat::rVector3d(0, a0(2), -a0(1) )  ): 
               atat::rVector3d( -a0(2) -a0(1), a0(0), a0(0) );

        // Then, lets create another... 
        a2( 0 ) = a0(1) * a1(2) - a0(2) * a1(1);
        a2( 1 ) = a0(2) * a1(0) - a0(0) * a1(2);
        a2( 2 ) = a0(0) * a1(1) - a0(1) * a1(0);
     }
    
     //! Strict weak ordering operator.
     bool LayerDepth :: operator()( const atat::rVector3d& _first, 
                                    const atat::rVector3d& _second ) const
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
     types::t_real LayerDepth :: operator()( const atat::rVector3d& _first ) const
     {
       __ASSERT( not isset, "Crystal::LayerDepth has not been set.\n" )
       const types::t_real a =  _first * a0; 
       return a - std::floor( a - types::tolerance );
     }

  } // namespace Crystal
} // namespace LaDa
