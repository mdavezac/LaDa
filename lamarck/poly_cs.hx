#ifndef _POLY_CS_H_
#define _POLY_CS_H_

#include <tinyxml/tinyxml.h>
#include <vector>

#include <opt/opt_polynome.h>

#include <lamarck/constituent_strain.h>
#include "poly_harmonic.h"

namespace VA_CE {

  #undef iVector3d 
  #define iVector3d atat::Vector3d<int>

  class Poly_CS : public opt::Polynome<double>
  {
    public:
      typedef opt::Polynome<double> :: TYPE TYPE;
      typedef opt::Polynome<double> :: CONTAINER CONTAINER;
      using opt::Polynome<double> :: evaluate_gradient;
      using opt::Polynome<double> :: evaluate_with_gradient;

    protected:
      using opt::Polynome<> :: variables;
      using opt::Polynome<> :: monomes;
      static const double ZERO_TOLERANCE;

    protected: 
      static std::vector<Poly_Harmonic> harmonics;
      
      // constructor, destructor, and helpers
    public:
      Poly_CS() : opt::Polynome<>() {};
      Poly_CS(const Poly_CS &_cs ) : opt::Polynome<>(_cs) {};
      virtual ~Poly_CS() {};
   
      bool operator = ( const Ising_CE::Constituent_Strain &_cs );

      bool Load_Harmonics( TiXmlElement *_el);
      virtual double evaluate() 
        { return opt::Polynome<> :: evaluate(); }
 
  };

} // namespace Ising_CE
#endif // _CONSTITTUENT_STRAIN_H_
