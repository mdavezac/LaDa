#include "LaDaConfig.h"

#include <algorithm>
#include <functional>

#include <opt/debug.h>
#include <math/misc.h>

#include "harmonic.h"

namespace LaDa
{
namespace CE
{
  namespace ConstituentStrain
  {
    namespace Harmonic
    {
      const std::string Cubic::type = "cubic";
      const std::string Tetragonal::type = "tetragonal";

      //! \cond
      namespace details
      {
        inline types::t_real pow2( types::t_real _t ) { return _t*_t; }
        inline types::t_real pow3( types::t_real _t ) { return _t * pow2( _t ); }
        inline types::t_real pow4( types::t_real _t ) { return pow2( pow2( _t ) ); }
        inline types::t_real pow6( types::t_real _t ) { return pow2( _t * pow2( _t ) ); }
        inline types::t_real pow8( types::t_real _t ) { return pow2( pow2( pow2( _t ) ) ); }
      }
      //! \endcond

      using math::ipow;
      // adds point in sorted list
      void Linear_Interpolator :: add_point (const types::t_real _x, const types::t_real _y )
      {
        if ( !points.size() )
        {
          points.push_back( Point(_x, _y) );
          return;
        }

        std::vector<Point> :: iterator i_point;
        i_point = std::find_if (points.begin(), points.end(),
                                std::bind2nd( std::mem_fun_ref( &Point::x_greater), _x ) );

        if ( i_point == points.end() )
          points.push_back( Point(_x, _y) );
        else if ( i_point->x == _x )
          return;
        else
          points.insert( i_point, Point(_x, _y) );
      }
      
      bool Linear_Interpolator :: Load (const TiXmlElement &_element)
      {
        double x=1.0, y=1.0;
        if ( not _element.Attribute("x", &x) )
          return false;
        if ( not _element.Attribute("y", &y) )
          return false;
        add_point(2.0*types::t_real(x)-1, types::t_real(y)); // on input, concentration is between 0 and 1
        return true;
      }

      types::t_real Linear_Interpolator :: evaluate( const types::t_real _x ) const
      {
        __ASSERT( points.size() == 0, "Nothing to interpolate\n" )

        std::vector<Point> :: const_iterator i_point;
        std::vector<Point> :: const_iterator i_begin = points.begin();
        std::vector<Point> :: const_iterator i_end = points.end();
        types::t_real p;

        i_point = std::find_if ( i_begin, i_end,
                                 std::bind2nd( std::mem_fun_ref( &Point::x_greater), _x ) );
     
        if ( i_point == i_begin )
          i_point++;
        if ( i_point == i_end )
          i_point--;
        i_begin = i_point - 1;
        p = ( i_point->y - i_begin->y ) / ( i_point->x - i_begin->x );
        return ( p * (_x - i_begin->x ) + i_begin->y ); 
      }
      types::t_real Linear_Interpolator :: evaluate_gradient( const types::t_real _x ) const
      {
        __ASSERT( points.size() == 0, "Nothing to interpolate\n" )

        std::vector<Point> :: const_iterator i_point;
        std::vector<Point> :: const_iterator i_begin = points.begin();
        std::vector<Point> :: const_iterator i_end = points.end();

        i_point = std::find_if ( i_begin, i_end,
                                 std::bind2nd( std::mem_fun_ref( &Point::x_greater), _x ) );
     
        if ( i_point == i_begin )
          i_point++;
        if ( i_point == i_end )
          i_point--;
        i_begin = i_point - 1;
        return ( ( i_point->y - i_begin->y ) / ( i_point->x - i_begin->x ) );
      }
      types::t_real Linear_Interpolator :: evaluate_with_gradient( const types::t_real _x,
                                                                   types::t_real &gradient ) const
      {
        __ASSERT( points.size() == 0, "Nothing to interpolate\n" )
        
        std::vector<Point> :: const_iterator i_point;
        std::vector<Point> :: const_iterator i_begin = points.begin();
        std::vector<Point> :: const_iterator i_end = points.end();

        i_point = std::find_if ( i_begin, i_end,
                                 std::bind2nd( std::mem_fun_ref( &Point::x_greater), _x ) );
     
        if ( i_point == i_begin )
          i_point++;
        if ( i_point == i_end )
          i_point--;
        i_begin = i_point - 1;
        gradient = ( i_point->y - i_begin->y ) / ( i_point->x - i_begin->x );
        return ( gradient * (_x - i_begin->x ) + i_begin->y ); 
      }


      types::t_real Cubic :: operator()( const math::rVector3d &_k ) const
      {
        types::t_real r2=_k.squaredNorm();
        types::t_real A=sqrt(1./(4.*M_PI));
        if ( math::leq( r2, types::t_real(0) ) )
        {
          switch ( rank )
          {
            case 0: return(0);
            case 1: return( -A*sqrt(21./4.)  );
            case 2: return( A*sqrt(13./2.)*(-5./4.) );
            case 3: return( 0.0 ); 
          }
        }
        switch (rank)
        {
           case 0: return(A);
           case 1: 
             return (   A*sqrt(21./4.)*(1.-5.*(math::ipow(_k(0),2)*ipow(_k(1),2) 
                      + ipow(_k(0),2)*ipow(_k(2),2) 
                      + ipow(_k(1),2)*ipow(_k(2),2)) / ipow(r2,2)));
           case 2: 
             return A * sqrt(13./2.) * (1./4.) * 
                    (   7 * (   ipow(_k(0),6) + ipow(_k(1),6) + ipow(_k(2),6)
                              + 30 * ipow(_k(0),2) * ipow(_k(1),2) * ipow(_k(2),2) 
                            ) / ipow(r2,3)
                      - 5
                    );
           case 3: 
             return A * sqrt(561.) * (1./8.) * 
                    (    ipow(_k(0),8) + ipow(_k(1),8) + ipow(_k(2),8) 
                       - 14. * (   ipow(_k(0),6)*ipow(_k(1),2)  
                                 + ipow(_k(0),6)*ipow(_k(2),2) 
                                 + ipow(_k(1),6)*ipow(_k(2),2) 
                                 + ipow(_k(0),2)*ipow(_k(1),6) 
                                 + ipow(_k(0),2)*ipow(_k(2),6) 
                                 + ipow(_k(1),2)*ipow(_k(2),6)
                               ) 
                       + 35. * (   ipow(_k(0),4)*ipow(_k(1),4) 
                                 + ipow(_k(0),4)*ipow(_k(2),4) 
                                 + ipow(_k(1),4)*ipow(_k(2),4)
                               )
                    ) / ipow(r2,4);
           default:
             __DOASSERT( true, "Kubic harmonics l > 3 not yet implemented.\n"  )
        }

        return 0.;
      }

      types::t_real Tetragonal :: operator()( const math::rVector3d &_k ) const
      {
        typedef types::t_real t_real;
        if( rank == 0 ) return t_real(1);
        types::t_real ctheta(1), cfai(1);
        types::t_real stheta(0), sfai(0);
        if(    math::neq(_k[0], t_real(0) )
            or math::neq(_k[1], t_real(0) ) )
        {
          ctheta =  details::pow2(_k[2]) / _k.squaredNorm();
          stheta =  t_real(1) - ctheta;
          types::t_real a = _k[0] * _k[0] + _k[1] * _k[1]; 
          cfai   = details::pow2( _k[0] ) / a;
          sfai   = details::pow2( _k[1] ) / a; 
        }
        if( rank == 1 ) return t_real(1.5) * ctheta - t_real(0.5);
        if( rank == 2 )
          return   (   t_real(35) * details::pow2( ctheta )
                     - t_real(30) * ctheta + t_real(3) ) 
                 * t_real( 0.125 );
        if( rank == 3 )
          return   ( details::pow2(cfai - sfai) - t_real(4) * cfai * sfai ) 
                 * details::pow2( stheta ) ;
        if( rank == 4 )
          return ( 
                    t_real(231) * details::pow3(ctheta)
                  - t_real(315) * details::pow2(ctheta) 
                  + t_real(105) * ctheta
                  - t_real(5)
                 ) * t_real( 0.0625 );
        if( rank == 5 )
        {
          return   ( t_real(11) * ctheta - t_real(1) )
                 * ( details::pow2(cfai - sfai) - t_real(4) * cfai * sfai )
                 * details::pow2( stheta ) 
                 * t_real(0.81675);
        }
        if( rank == 6 )
          return (   
                     t_real( 6435) * details::pow4( ctheta )
                   - t_real(12012) * details::pow3( ctheta ) 
                   + t_real( 6930) * details::pow2( ctheta ) 
                   - t_real( 1260) * ctheta 
                   + t_real(   35) 
                 ) * std::sqrt( t_real(17) ) / t_real( 128 );
        if( rank == 7 )
          return   ( details::pow2(cfai - sfai) - t_real(4) * cfai * sfai )
                 * details::pow2(stheta) 
                 * t_real(3) / t_real(64) * std::sqrt( t_real( 654.5 ) )
                 * (   t_real(65) * details::pow2(ctheta) 
                     - t_real(26) * ctheta
                     + t_real(1) );
        if( rank == 8 )
        {
          t_real c = t_real(1) - t_real(8) * cfai * sfai; c *= c; 
          t_real d = cfai - sfai; d *= d * t_real(8) * cfai * sfai;
          return   ( c - d ) // cos 8fai
                 * details::pow4( stheta ) 
                 * t_real(3) / t_real(128) * std::sqrt( t_real(6077.5) );
        }
                   

        __DOASSERT(true, "rank > 8 not implemented for tetragonal harmonics.\n" )
        return t_real(0);
      }

      
    } // namespace Harmonic
  } // namespace ConstituentStrain

} // namespace CE
}
