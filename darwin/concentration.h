//
//  Version: $Id$
//
#ifndef _CONCENETRATION_H_
#define _CONCENETRATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdexcept>       // std::runtime_error
#include "lamarck/structure.h"
#include "opt/types.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif


  class X_vs_y
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<X_vs_y>( X_vs_y & );
#endif
    protected:
      types::t_real a, b, c;
      types::t_real x0;
      types::t_real y0;
      bool singlec;

    public:
      X_vs_y() : a(0), b(0), c(0), x0(0), y0(0), singlec(false) {}
      X_vs_y( const X_vs_y &_c) : a(_c.a), b(_c.b), c(_c.c), 
                                  x0(_c.x0), y0(_c.y0), singlec(_c.singlec) {}
      bool Load( const TiXmlElement &_node );
      types::t_real get_x( types::t_real _y ) 
      {
        if ( singlec ) return x0;
        return c + b * _y + a * _y * _y; 
      }
      void set_xy( types::t_real _x, types::t_real _y )
      {
        x0 = _x; y0 = _y; singlec = true;
      } 
      types::t_real get_y() { return y0; }
      types::t_real get_x() { return x0; }
      types::t_real get_y( types::t_real _x )
      {
        if ( singlec ) return y0;
        if ( std::abs ( a ) < types::tolerance )
          return ( _x - c ) / b;
       
        types::t_real det = b*b - 4.0 * (c-_x) * a; 
        if ( det < 0 )
        {
          std::cerr << "Error when using Concentration::get_y(" << _x<<")" << std::endl
                    << "determinent is negative, " << det << std::endl;
          throw std::runtime_error("");
        }
        det = std::sqrt(det);
        types::t_real u = 1.0 / ( 2.0 * a );  
        types::t_real r0 =  (-b + det ) * u;
        types::t_real r1 =  (-b - det ) * u;
        if ( std::abs(r0 - 1.0 ) < types::tolerance ) r0 = 1.0;
        if ( std::abs(r0 + 1.0 ) < types::tolerance ) r0 = -1.0;
        if ( std::abs(r1 - 1.0 ) < types::tolerance ) r1 = 1.0;
        if ( std::abs(r1 + 1.0 ) < types::tolerance ) r1 = -1.0;
        if (     ( r0 < -1.0 or r0 > 1.0 )
             and ( r1 < -1.0 or r1 > 1.0 ) )
        {
          std::cerr << a + b +c << " " << a - b + c << std::endl;
          std::cerr << "Error when using Concentration::get_y(" << _x<< ")" << std::endl;
          std::cerr << " r0= " << r0  << " and r1= " << r1 << std::endl;
          throw std::runtime_error("");
        }
        if ( r0 < -1.0 or r0 > 1.0 ) 
          return r1; 
        return r0;
      }
      bool can_inverse( types::t_real _x );
      bool is_singlec () const { return singlec; }

  };

types::t_real set_concentration( Ising_CE::Structure &_str,
                                 types::t_real _target = -2.0);


#endif 
