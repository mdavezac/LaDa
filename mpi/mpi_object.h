//
//  Version: $Id$
//
#ifndef _MPI_OBJECT_H_
#define _MPI_OBJECT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "macros.h"

#ifdef _MPI

#ifdef _MPICH_MPI_
#include <mpi2c++/mpi++.h>
#elif defined(_OPENMPI_MPI_)
#include <mpi.h>
#endif

#include <stdexcept>       // std::runtime_error
#include <iostream>
#include <math.h>

#include <opt/types.h>

#include "base.h"
#include "comm.h"
#include "standalone.h"


namespace mpi
{
  //! \brief Main object for creation and destruction of mpich, as well as for
  //!        intialization of transaction helpers
  extern InitDestroy main;
}

#endif
#endif
