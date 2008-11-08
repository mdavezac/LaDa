//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "mpi_object.h"

#ifdef _MPI
namespace LaDa
{
  namespace mpi
  {
    boost::mpi::communicator *main = NULL; 
  }
} // namespace LaDa
#endif
