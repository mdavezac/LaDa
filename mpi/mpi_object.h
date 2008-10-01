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
#include <complex>

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
    //! Serializes a complex number.
    template<class Archive, class T>
    inline void serialize (Archive &ar, std::complex<T>& z, const unsigned int file_version)
    {
      ar & boost::serialization::make_nvp ("real", real(z));
      ar & boost::serialization::make_nvp ("imag", imag(z));
    }

  }

}

namespace mpi
{
  //! \brief Main object for creation and destruction of mpich, as well as for
  //!        intialization of transaction helpers
  extern boost::mpi::communicator *main;

  //! Adds a communicator pointer to a class.
  class AddCommunicator
  {
    public:
      //! Constructor.
      AddCommunicator() : comm_(NULL) {} 
      //! Constructor.
      AddCommunicator( boost::mpi::communicator &_c ) : comm_( &_c ) {} 
      //! Copy Constructor.
      AddCommunicator( const AddCommunicator &_c ) : comm_( _c.comm_ ) {} 

      //! Sets mpi pointer.
      void set_mpi( boost::mpi::communicator* _c )
      {
        __ASSERT( _c == NULL, "Pointer not set.\n" )
        comm_ = _c;
      }
      
      //! Returns reference to communicator.
      boost::mpi::communicator &comm()
      {
        __ASSERT( comm_ == NULL, "Pointer not set.\n" )
        return *comm_;
      }
      //! Returns a constant reference to communicator.
      const boost::mpi::communicator &comm() const 
      {
        __ASSERT( comm_ == NULL, "Pointer not set.\n" )
        return *comm_;
      }
    protected:

      //! The MPI Communicator.
      boost::mpi::communicator *comm_;
  };
}

#define MPI_COMMDEC ::mpi::AddCommunicator
#define MPI_COMMA ,
#define MPI_COMMCOPY( _c ) ::mpi::AddCommunicator( _c )
#define MPI_COMM ::mpi::AddCommunicator::comm()
#define MPI_FORWARD_MEMBERS( base ) \
  //! Allows derived classes to have access to ::mpi::AddCommunicator members. \
  boost::mpi::communicator &comm() { return base::comm(); } \
  //! Allows derived classes to have access to ::mpi::AddCommunicator members. \
  const boost::mpi::communicator &comm() const { return base::comm(); } \
  //! Allows derived classes to have access to ::mpi::AddCommunicator members. \
  void set_mpi( boost::mpi::communicator* _c ) { base::set_mpi( c ); } \

#else

#define MPI_COMMDEC
#define MPI_COMMA
#define MPI_COMMCOPY( _c ) 
#define MPI_GETCOMM
#define MPI_FORWARD_MEMBERS( base ) 

#endif

#endif
