#ifndef _MPI_MACROS_H_
# define _MPI_MACROS_H_

#include "LaDaConfig.h"

# ifndef LADA_MPI
#   define LADA_MPI_CODE(code) 
#   define LADA_ROOT(comm, code) code
#   define LADA_MPI_START
    
# else

#   include <opt/debug.h>
#   define LADA_MPI_CODE(code) code
#   define LADA_ROOT(__comm, code) if( __comm.rank() == 0 ) { code } 
#   define LADA_MPI_START \
         boost::scoped_ptr< boost::mpi::environment > env; \
         boost::scoped_ptr< boost::mpi::communicator > world; \
         try \
         { \
           env.reset( new boost::mpi::environment (argc, argv) ); \
           world.reset( new boost::mpi::communicator ); \
         } \
         catch ( std::exception &_e ) \
         { \
           std::cerr << "Error encountered while creating MPI environment.\n"; \
           return 0; \
         } \

# endif
#endif
