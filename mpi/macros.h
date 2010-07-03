#ifndef _MPI_MACROS_H_
# define _MPI_MACROS_H_

#include "LaDaConfig.h"

# ifndef _MPI
#   define __MPICODE(code) 
#   define __MPICONSTRUCTORCODE(code) 
#   define __TRYMPICODE(code, error) 
#   define __MPIROOT(comm, code) 
#   define __ROOTCODE(comm, code) code
#   define __NOTMPIROOT(comm, code) 
#   define __SERIALCODE(code) code
#   define __MPISEQUENTIAL(code) 
#   define __DOMPISEQUENTIAL(code) code
#   define __DOCOMMSEQUENTIAL(code, comm) code
#   define __MPISERIALCODE(coda, codb) codb
#   define __MPI_START__
    
    // mpi::Broadcast macros
#   define  ___DECLAREMPIOBJECT( name )
#   define ___FRIENDMPIOBJECT( name )


# else

#   include <opt/debug.h>
#   define __MPICODE(code) code
#   define __MPICONSTRUCTORCODE(code) , code
#   define __TRYMPICODE(code, error) try { code }\
           catch ( std::exception &_e )\
           {\
             std::ostringstream sstr;\
             sstr << __SPOT_ERROR << error << _e.what();\
             throw std::runtime_error( sstr.str() );\
           }
#   define __MPIROOT(__comm, code) if( __comm.rank() == 0 ) { code }
#   define __ROOTCODE(__comm, code) __MPIROOT( __comm, code)
#   define __NOTMPIROOT(__comm, code) if( __comm.rank() != 0 ) { code }
#   define __SERIALCODE(code) 
#   define __DOCOMMSEQUENTIAL(comm, code) \
       for( types::t_int i = 0; i < comm.size(); ++i )\
       {\
         if ( comm.rank() == i ) { code }\
           comm.barrier();\
       }
#   define __DOMPISEQUENTIAL(comm, code)  __DOCOMMSEQUENTIAL( comm, code )
#   define __MPISEQUENTIAL(code) __DOMPISEQUENTIAL(code)
#   define __MPISERIALCODE(coda, codb) coda
#   define __MPI_START__ \
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
