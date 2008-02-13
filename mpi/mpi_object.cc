//
//  Version: $Id$
//
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <stdexcept>       // std::runtime_error

#include <atat/vectmac.h>
#include <opt/debug.h>

#include "mpi_object.h"

#ifdef _MPI
namespace mpi
{
  InitDestroy main; 
  const types::t_int ROOT_NODE = 0;
}

#endif
