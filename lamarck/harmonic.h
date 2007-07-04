#ifndef _HARMONICS_H_
#define _HARMONICS_H_


#include <vector>
#include <math.h>

#include <tinyxml/tinyxml.h>


#include "opt/types.h"
#include "atat/vectmac.h"

#ifdef _MPI 
  #include "mpi/mpi_object.h"
#endif

namespace Ising_CE 
{

  class Linear_Interpolator 
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<Linear_Interpolator> ( Linear_Interpolator& );
#endif
    static const types::t_real tolerance;
    
    struct Point
    {
      types::t_real x,y;
      Point() {}
      Point(const Point &_p) { x = _p.x; y = _p.y; }
      Point(const types::t_real _x, const types::t_real _y) { x = _x; y = _y; }
      void operator= (const Point &_p) { x = _p.x; y = _p.y; }
      bool x_greater( const types::t_real _x ) const
        { return (x > _x); }
      bool x_lesser( const types::t_real _x ) const
        { return (x < _x); }
    };

    std::vector<Point> points;

    public:
      Linear_Interpolator(){};
      Linear_Interpolator(const Linear_Interpolator &_inter)
        { points = _inter.points; };

      void operator=(const Linear_Interpolator &_inter)
        { points = _inter.points; }; 

      void add_point( const types::t_real _x, const types::t_real _y);
      types::t_real evaluate( const types::t_real _x ) const;
      types::t_real evaluate_gradient( const types::t_real _x ) const;
      types::t_real evaluate_with_gradient( const types::t_real _x, types::t_real &gradient ) const;
      bool Load (const TiXmlElement &_element);
      void clear()
        { points.clear(); };
  };


  class Harmonic 
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<Harmonic> ( Harmonic& );
#endif
    static const types::t_real tolerance;

    Linear_Interpolator interpolation;
    types::t_unsigned rank;
    static types::t_real attenuation;

    public:
      Harmonic() {};
      Harmonic( const Harmonic &_h )
        { rank = _h.rank; interpolation = _h.interpolation; }

      void operator=(const Harmonic &_h )
        { rank = _h.rank; interpolation = _h.interpolation; }

      // behavior required by minimizer::??
    public:
      inline void set_variable(const types::t_unsigned i, const types::t_real term);
      inline types::t_real get_variable(const types::t_unsigned i) const;
      inline types::t_unsigned get_nb_variables() const;
      inline types::t_real evaluate(const types::t_real _x, const atat::rVector3d &vec) const;
      inline types::t_real evaluate(const types::t_real _x) const;
      inline types::t_real evaluate(const atat::rVector3d &vec) const;
      inline types::t_real evaluate_gradient(const types::t_real _x, const atat::rVector3d &vec) const;
      inline types::t_real evaluate_gradient(const types::t_real _x) const;
      inline types::t_real evaluate_with_gradient(const types::t_real _x,
                                                  const atat::rVector3d &vec, 
                                                  types::t_real &gradient) const;
      inline types::t_real evaluate_with_gradient(const types::t_real _x,
                                                  types::t_real &gradient) const;

      // other
    public:
      types::t_real evaluate_harmonic(const atat::rVector3d &k) const;
      
      bool Load(const TiXmlElement &_element);
      static void set_attenuation( const types::t_real _a ) 
        { attenuation = (_a == 0) ? 0 : 1.0 / (_a*_a); }
      void clear()
        { interpolation.clear(); };

  };

  inline types::t_real Harmonic :: evaluate(const types::t_real _x, const atat::rVector3d &k) const
  {
    return (   interpolation.evaluate(_x) 
             * exp( -norm2(k) * attenuation )
             * evaluate_harmonic( k ) );
  }
  inline types::t_real Harmonic :: evaluate(const types::t_real _x) const
  {
    return (  interpolation.evaluate(_x) );
  }
  inline types::t_real Harmonic :: evaluate(const atat::rVector3d &k) const
  {
    return (  exp( -norm2(k) * attenuation )
             * evaluate_harmonic( k ) );
  }
  inline types::t_real Harmonic :: evaluate_gradient(const types::t_real _x, const atat::rVector3d &k) const
  {
    return (   interpolation.evaluate_gradient(_x) 
             * exp( -norm2(k) *  attenuation )
             * evaluate_harmonic( k ) );
  }
  inline types::t_real Harmonic :: evaluate_gradient(const types::t_real _x) const
  {
    return ( interpolation.evaluate_gradient(_x) );
  }
  inline types::t_real Harmonic :: evaluate_with_gradient(const types::t_real _x, 
                                                   const atat::rVector3d &k,
                                                   types::t_real &gradient) const
  {
    types::t_real factor = exp( -norm2(k) * attenuation )
                    * evaluate_harmonic( k );
    types::t_real result =   interpolation.evaluate_with_gradient(_x, gradient) 
                    * factor ;
    gradient *= factor;
    return result;
  }
  inline types::t_real Harmonic :: evaluate_with_gradient(const types::t_real _x,
                                                   types::t_real &gradient) const
  {
    return interpolation.evaluate_with_gradient(_x, gradient);
  }

} // namespace Ising_CE 

#endif // _HARMONICS_H_
