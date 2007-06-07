#ifndef _POLY_HARMONICS_H_
#define _POLY_HARMONICS_H_

#include <vector>
#include <utility>
#include <utility>
#include <math.h>
#include <tinyxml/tinyxml.h>
#include <atat/vectmac.h>
#include <opt/opt_polynome.h>



namespace VA_CE 
{

  #undef rVector3d 
  #define rVector3d atat::Vector3d<Real>
  class Poly_Harmonic
  {
    static const double ZERO_TOLERANCE;

    std::vector< std::pair<double, int> > poly_coefs;
    unsigned rank;
    static double attenuation;

    public:
      Poly_Harmonic() {};
      Poly_Harmonic( const Poly_Harmonic &_h )
        { rank = _h.rank; poly_coefs = _h.poly_coefs; }

      void operator=(const Poly_Harmonic &_h )
        { rank = _h.rank; poly_coefs = _h.poly_coefs; }

      // behavior required by opt::Minimize
    public:
      inline double evaluate(const double _x, const rVector3d &vec) const;
      double evaluate(const double _x) const;
      inline double evaluate(const rVector3d &vec) const;
      inline double evaluate_gradient(const double _x, const rVector3d &vec) const;
      double evaluate_gradient(const double _x) const;
      inline double evaluate_with_gradient(const double _x, const rVector3d &vec, double &gradient) const;
      double evaluate_with_gradient(const double _x, double _grad) const;

      // other
    public:
      double evaluate_harmonic(const rVector3d &k) const;
      
      bool Load( TiXmlElement *element);
      static void set_attenuation( const double _a ) 
      { 
        _a == 0 ?  attenuation = 0 : attenuation = 4.0 / (_a*_a);
      }
      void clear()
        { poly_coefs.clear(); };

      void var_change( opt::Polynome<> &_poly, bool _linearize = false ) const;
  };

  inline double Poly_Harmonic :: evaluate(const double _x, const rVector3d &k) const
  {
    return (   evaluate(_x) 
             * exp( -norm2(k) * attenuation )
             * evaluate_harmonic( k ) );
  }
  inline double Poly_Harmonic :: evaluate(const rVector3d &k) const
  {
    return (  exp( -norm2(k) * attenuation )
             * evaluate_harmonic( k ) );
  }
  inline double Poly_Harmonic :: evaluate_gradient(const double _x, const rVector3d &k) const
  {
    return (   evaluate_gradient(_x) 
             * exp( -norm2(k) * attenuation )
             * evaluate_harmonic( k ) );
  }
  inline double Poly_Harmonic :: evaluate_with_gradient(const double _x, 
                                                   const rVector3d &k,
                                                   double &gradient) const
  {
    double factor = exp( -norm2(k) * attenuation )
                    * evaluate_harmonic( k );
    double result =   evaluate_with_gradient(_x, gradient) 
                    * factor ;
    gradient *= factor;
    return result;
  }

} // namespace Ising_CE 

#endif // _HARMONICS_H_
