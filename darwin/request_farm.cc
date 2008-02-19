//
//  Version: $Id$
//

#ifdef _MPI

#include "request_farm.h"


namespace GA
{
  namespace mpi
  {
    Bull :: Bull   ( MPI::Comm &_farmercomm, MPI::Comm & _cowcomm)
                 : ::mpi::Base( _farmercomm ), cowcomm(NULL),
                    to_cows( NULL ), farmer_in(t_FromFarmer::UNDEFINED),
                    farmer_out(t_ToFarmer::UNDEFINED), cows_out(NULL), ncows(0)
   {
      __ASSERT( not ::mpi::main::is_root_node() ); 
      __ASSERT( (not _cowcomm) or _cowcomm->rank() == ::mpi::ROOT_NODE ); 
   }

}


#endif
