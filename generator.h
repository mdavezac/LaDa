#ifndef _GENERATOR_H_
#define _GENERATOR_H_

#include <eo/utils/eoRndGenerators.h>

namespace LaDa 
{
  class Generator : public eoRndGenerator<double>
  {
    private: 
      eoRng& uniform;
      
    public:
      Generator( eoRng &_rng = rng) : uniform(_rng) {};

      double operator()(void)
        { return ( uniform.flip() ) ? -1.0 : 1.0; }
  };

} // namespace LaDa

#endif // _GENERATOR_H_
