//
//  Version: $Id$
//
#ifndef _MPI_MACROS_H_
#define _MPI_MACROS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef _MPI
#define __MPICODE(code) 
#define __MPICONSTRUCTORCODE(code) 
#define __TRYMPICODE(code, error) 
#define __MPIROOT(code) 
#define __ROOTCODE(code) code
#define __NOTMPIROOT(code) 
#define __SERIALCODE(code) code
#define __MPISEQUENTIAL(code) 
#define __DOMPISEQUENTIAL(code) code
#define __DOCOMMSEQUENTIAL(code, comm) code
#define __MPISERIALCODE(coda, codb) codb

// mpi::Broadcast macros
#define  ___DECLAREMPIOBJECT( name )
#define ___FRIENDMPIOBJECT( name )
#else
#include <opt/debug.h>
#define __MPICODE(code) code
#define __MPICONSTRUCTORCODE(code) , code
#define __TRYMPICODE(code, error) try { code }\
        catch ( std::exception &_e )\
        {\
          std::ostringstream sstr;\
          sstr << __SPOT_ERROR << error << _e.what();\
          throw std::runtime_error( sstr.str() );\
        }
#define __MPIROOT(code) if( ::mpi::main.rank() != 0 ) { code }
#define __ROOTCODE(code) __MPIROOT(code)
#define __NOTMPIROOT(code) if( not ::mpi::main.rank() != 0 ) { code }
#define __SERIALCODE(code) 
#define __DOCOMMSEQUENTIAL(comm, code) \
    for( types::t_int i = 0; i < comm.size(); ++i )\
    {\
      if ( comm.rank() == i ) { code }\
        comm.barrier();\
    }
#define __DOMPISEQUENTIAL(code)  __DOCOMMSEQUENTIAL( ::mpi::main, code )
#define __MPISEQUENTIAL(code) __DOMPISEQUENTIAL(code)
#define __MPISERIALCODE(coda, codb) coda

// mpi::Broadcast macros
#define ___DECLAREMPIOBJECT( name )\
 namespace mpi \
 { \
   template<> class BroadCast::Object< name >;\
 }
#define ___FRIENDMPIOBJECT( name )\
  friend class mpi::BroadCast::Object< name >;
#endif
#endif
