//
//  Version: $Id$
//
#ifndef _DARWIN_GENCOUNT_H_
#define _DARWIN_GENCOUNT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "opt/types.h"

namespace darwin
{
  // generation counter
  class GenCount
  {
    protected:
      types::t_unsigned age;
    public:
      GenCount( const GenCount &_gc) : age(_gc.age) {};
      GenCount( types::t_unsigned _age) : age(_age) {};
      GenCount() : age(0) {};
      void operator ++() 
        { ++age; }
      types::t_unsigned operator()() const
        { return age; }
  };
} // namespace LaDa

#endif
