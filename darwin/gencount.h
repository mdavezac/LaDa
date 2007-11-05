//
//  Version: $Id$
//
#ifndef _DARWIN_GENCOUNT_H_
#define _DARWIN_GENCOUNT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "opt/types.h"

namespace GA
{
  //! \brief generation counter.
  //! \details This is simply something which can be incremented and only
  //!          incremented by one. It should be so once every generation. It
  //!          can also return its current value.
  class GenCount
  {
    protected:
      //! the number of generation since start of %GA run
      types::t_unsigned age;
    public:
      //! Copy Constructor
      GenCount( const GenCount &_gc) : age(_gc.age) {};
      //! Constructor and initializer
      GenCount( types::t_unsigned _age) : age(_age) {};
      //! Constructor
      GenCount() : age(0) {};
      //! increments the number of generations by one
      void operator ++() 
        { ++age; }
      //! Returns the number of generations
      types::t_unsigned operator()() const
        { return age; }
  };
} // namespace GA

#endif
