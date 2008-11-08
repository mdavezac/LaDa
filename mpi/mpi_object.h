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
#include <boost/serialization/complex.hpp>
#include <boost/type_traits/add_reference.hpp>
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
    void serialize(Archive & ar, LaDa::atat::rVector3d & g, const unsigned int version)
     { ar & g.x; }
     //! Serializes atat integer vectors.
    template<class Archive>
    void serialize(Archive & ar, LaDa::atat::iVector3d & g, const unsigned int version)
     { ar & g.x; }
    //! Serializes atat real matrices.
    template<class Archive>
    void serialize(Archive & ar, LaDa::atat::rMatrix3d & g, const unsigned int version)
     { ar & g.x; }
    //! Serializes atat integer matrices.
    template<class Archive>
    void serialize(Archive & ar, LaDa::atat::iMatrix3d & g, const unsigned int version)
     { ar & g.x; }
  }

}

namespace LaDa
{
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

    //! A class to easily create and pass two object to mpi.
    template< class FIRST, class SECOND >
    struct Pair
    {
      //! A reference type to the first type.
      typedef typename boost::add_reference<FIRST> :: type t_First;
      //! A reference type to the second type.
      typedef typename boost::add_reference<SECOND> :: type t_Second;
      //! The first reference.
      t_First first;
      //! The second reference.
      t_Second second;
      //! The constructor and initializer.
      Pair( t_First _first, t_Second _second ) : first( _first ), second( _second ) {}
      //! The copy constructor.
      Pair( const Pair &_c ) : first( _c.first ), second( _c.second ) {}
      //! The serialize member.
      template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version)
        { _ar & first; _ar & second; }
    };
  }
}

#define MPI_COMMDEC ::LaDa::mpi::AddCommunicator
#define MPI_COMMA ,
#define MPI_COMMCOPY( _c ) ::LaDa::mpi::AddCommunicator( _c )
#define MPI_COMM ::LaDa::mpi::AddCommunicator::comm()
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
