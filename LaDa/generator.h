#ifndef _GENERATOR_H_
#define _GENERATOR_H_

#include <eo/utils/eoRndGenerators.h>
#include <opt/types.h>

namespace LaDa 
{
  class Generator : public eoRndGenerator<types::t_real>
  {
    private: 
      eoRng& uniform;
      
    public:
      Generator( eoRng &_rng = rng) : uniform(_rng) {};

      types::t_real operator()(void)
        { return ( uniform.flip() ) ? -1.0 : 1.0; }
  };

} // namespace LaDa

#endif // _GENERATOR_H_
