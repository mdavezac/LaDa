#include "LaDaConfig.h"

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
