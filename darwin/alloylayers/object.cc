//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <algorithm>

#include <opt/random.h>
#include <opt/types.h>

#include "object.h"

namespace LaDa
{
  namespace GA
  {
    namespace AlloyLayers 
    {
      struct comp
      {
        bool operator()( types::t_real const &_a, types::t_real const &_b ) const
          { return std::abs(_a-_b) < 1e-12; }
      };
      bool Object :: operator==( Object const &_c ) const
      {
        for(size_t start(0), end(_c.bitstring.size()); start < end; ++start)
        {
          if( std::equal(bitstring.begin()+start, bitstring.end(),
                         _c.bitstring.begin()+start, comp()) )
          {
            if( std::equal(bitstring.begin(), bitstring.begin()+start,
                           _c.bitstring.begin(), comp()) ) return true;
          }
          if( std::equal(bitstring.rbegin()+start, bitstring.rend(),
                         _c.bitstring.rbegin()+start, comp()) )
          {
            if( std::equal(bitstring.rbegin(), bitstring.rbegin()+start,
                           _c.bitstring.rbegin(), comp()) ) return true;
          }
        }
        return false;
      }

      bool Object::random_symmetry()
      {
        t_Container copy( bitstring );
        size_t const size(bitstring.size());
        size_t const start(eo::rng.uniform(size));
        if(  eo::rng.flip() )
        {
          std::copy(copy.begin()+start, copy.end(), bitstring.begin());
          std::copy(copy.begin(), copy.begin()+start, bitstring.end()-start);
        }
        else
        {
          std::copy(copy.rbegin()+start, copy.rend(), bitstring.begin());
          std::copy(copy.rbegin(), copy.rbegin()+start, bitstring.end()-start);
        }
      }
    } // namespace Keepers
  } // namespace GA
}


