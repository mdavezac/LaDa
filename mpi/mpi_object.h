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
#include <boost/mpi/communicator.hpp>

// #ifdef _MPICH_MPI_
// #include <mpi2c++/mpi++.h>
// #elif defined(_OPENMPI_MPI_)
// #include <mpi.h>
// #endif

// #include <stdexcept>       // std::runtime_error
// #include <iostream>
// #include <math.h>

#include <atat/vectmac.h>
#include <opt/types.h>

// #include "base.h"
// #include "comm.h"
// #include "standalone.h"

namespace boost {
  namespace serialization {

    //! Serializes atat real vectors.
    template<class Archive>
    void serialize(Archive & ar, atat::rVector3d & g, const unsigned int version)
     { ar & g.x; }
     //! Serializes atat integer vectors.
    template<class Archive>
    void serialize(Archive & ar, atat::iVector3d & g, const unsigned int version)
     { ar & g.x; }
    //! Serializes atat real matrices.
    template<class Archive>
    void serialize(Archive & ar, atat::rMatrix3d & g, const unsigned int version)
     { ar & g.x; }
    //! Serializes atat integer matrices.
    template<class Archive>
    void serialize(Archive & ar, atat::iMatrix3d & g, const unsigned int version)
     { ar & g.x; }

  }

}

namespace mpi
{
  //! \brief Main object for creation and destruction of mpich, as well as for
  //!        intialization of transaction helpers
  extern boost::mpi::communicator *main;
}

#endif
#endif
