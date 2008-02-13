//
//  Version: $Id$
//
#ifndef _MPI_MACROS_H_
#define _MPI_MACROS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef _MPI
#define __DOMPICODE(code) 
#define __TRYMPICODE(code, error) 
#define __MPIROOT(code) 
#define __ROOTCODE(code) code
#define __NOTMPIROOT(code) 
#define __SERIALCODE(code) code
#define __MPISEQUENTIAL(code) 
#define __MPISERIALCODE(coda, codb) codb
#else
#include <opt/debug.h>
#define __DOMPICODE(code) code
#define __TRYMPICODE(code, error) try { code }\
        catch ( std::exception &_e )\
        {\
          std::ostringstream sstr;\
          sstr << __SPOT_ERROR << error << _e.what();\
          throw std::runtime_error( sstr.str() );\
        }
#define __MPIROOT(code) if( mpi::main.is_root_node() ) { code }
#define __ROOTCODE(code) __MPIROOT(code)
#define __NOTMPIROOT(code) if( not mpi::main.is_root_node() ) { code }
#define __SERIALCODE(code) 
#define __MPISEQUENTIAL(code) \
    for( types::t_int i = 0; i < mpi::main.size(); ++i )\
    {\
      if ( mpi::main.rank() == i ) { code }\
      mpi::main.barrier();\
    }
#define __MPISERIALCODE(coda, codb) coda
#endif
