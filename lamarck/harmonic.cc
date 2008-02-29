//
//  Version: $Id$
//
#include <algorithm>
#include <functional>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/debug.h>

#include "harmonic.h"

namespace Ising_CE
{
  types::t_real Harmonic :: attenuation = 1000000.0;
  using atat::ipow;
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
  types::t_real Linear_Interpolator :: evaluate_with_gradient( const types::t_real _x, types::t_real &gradient ) const
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


  types::t_real Harmonic :: evaluate_harmonic ( const atat::rVector3d &_k ) const
  {
    types::t_real r2=norm2(_k);
    types::t_real A=sqrt(1./(4.*M_PI));
    if ( fabs(r2) <= types::tolerance )
    {
      switch ( rank )
      {
        case 0:
          return(0);
        case 1:
          return( -A*sqrt(21./4.)  );
        case 2:
          return( A*sqrt(13./2.)*(-5./4.) );
        case 3:
          return( 0.0 ); 
      }
    }
    switch (rank)
    {
       case 0:
         return(A);
       case 1: 
         return (   A*sqrt(21./4.)*(1.-5.*(ipow(_k(0),2)*ipow(_k(1),2) 
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
         std::cerr << "Kubic harmonics l > 3 not yet implemented"
                   << std::endl;
         exit(0);
    }

    return 0.;
  }
  
  bool Harmonic :: Load (const TiXmlElement &_element) 
  {
    const TiXmlElement *child;
    int i=1;

    
    if ( not _element.Attribute("rank", &i) )
    {
      std::cerr << "Harmonic has no rank on input"
                << std::endl;
      return false;
    }
    if ( i < 0 || i > 3)
    {
      std::cerr << "rank of harmonic is out of range on input"
                << std::endl;
      return false;
    }
    rank = types::t_unsigned(abs(i));
    child = _element.FirstChildElement( "Point" );
    clear(); // clear interpolation
    if ( !child )
      return false;
    for ( ; child; child = child->NextSiblingElement( "Point" ) )
      if ( not interpolation.Load( *child ) )
        return false;

    return true;
  }

} // namespace Ising_CE

