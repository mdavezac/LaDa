//
//  Version: $Id$
//
#ifndef _MPI_COMM_H_
#define _MPI_COMM_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <opt/types.h>
#include "base.h"

namespace mpi
{
  //! \brief %Base class for all data communication classes
  //! \details Defines integer, real, and char buffers to be used for data transactions.
  //! It cannot be used as is, but should be derived from.
  class CommBase : public Base
  {
    protected:
      types::t_int *int_buff,       //!< beginning of integer buffer
                   *end_int_buff,   //!< end of integer buffer
                   *cur_int_buff;   //!< current position in the integer buffer
      types::t_char *char_buff,     //!< beginning of character buffer
                    *end_char_buff, //!< end of character buffer
                    *cur_char_buff; //!< current position in the character buffer
      types::t_real *real_buff,     //!< beginning of real buffer
                    *end_real_buff, //!< end of real buffer                  
                    *cur_real_buff; //!< current position in the real buffer
      //! \brief array to communicate size of all three buffers
      types::t_unsigned buffer_size[3]; 


    public:
      //! Constructor, initializes all variables to NULL or zero
      CommBase ();
      //! \brief initialization Constructor, copies only mpi::Base members.
      //! \details Should be used in conjunction with mpi::main to initialize any data
      //! communication helper class. Only mpi::Base members are copied,
      //! Members from CommBase are set to 0 or NULL. 
      CommBase   ( const Base &_c );
      //! \brief true Copy Constructor
      CommBase   ( const CommBase &_c);
      //! \brief Destructor. Does nothing ;)
      //! \details Destructor leaves it to derived classes to destroy buffers!!
      virtual ~CommBase() {}

      //! Destroys all buffers if they exist, then assigns all members to 0 or NULL;
      void destroy_buffers();
      //! \brief Allocates buffers
      //! \details should be overidden in
      //! derived classes in order to create buffers according to data
      //! transaction scheme and CommBase::buffer_size. May become non-virtual
      //! in the future.
      //! \param _root optionally, specifies root process of data transaction
      virtual bool allocate_buffers( types::t_unsigned _root = ROOT_NODE ) = 0;
      //!\brief  Historycally, was meant to be the main function of all data
      //! transaction helper classes. 
      //! \details Historycally, was meant to be the main interface function of all
      //! data transaction helper classes. It is insufficient in itself (see
      //! BroadCast::serialize() below), and should be overidden with
      //! mpi::operator<<() when possible
      //! \param _root optionally, specifies root process of data transaction
      virtual bool operator()( types::t_unsigned _root = ROOT_NODE ) = 0;
  };

  //! \brief Can be used to send data from one root process to all processes.
  //! \details It adopts a six fold scheme:
  //!  - all variables should be sized up using BroadCast::serialize(), in
  //!      order to determine the size of the buffers. Any number of calls to 
  //!      BroadCast::serialize() can be made, with any type of argument for
  //!      which it has been explicitely defined,
  //!  - then the buffers should be allocated using BroadCast::allocate_buffers(),
  //!  - a second series of calls to BroadCast::serialize should be made in
  //!      order to load up the buffers. 
  //!  - the data in the buffers should be broadcasted using Broadcast::operator()(),
  //!  - the data should be forwarded from the buffers to the serizalized
  //!      object using BroadCast::serialize(), in exactly the same order as
  //!      they were loaded in.
  //!  - a last optional call to BroadCast::operator()() readies the object for
  //!      another data transaction.
  //!  .
  //! For instance, we have
  //! \code
  //!  mpi::BroadCast bc(mpi::main); // initializes a BroadCast objet named bc
  //!  bc.serialize(A);  // sizes up object A
  //!  bc.serialize(B);  // sizes up object B
  //!  bc.serialize(C);  // sizes up object C
  //!  bc.allocate_buffers(); // allocates buffers
  //!  bc.serialize(A);  // loads object A into buffers
  //!  bc.serialize(B);  // loads object B into buffers
  //!  bc.serialize(C);  // loads object C into buffers
  //!  bc(); // broadcasts buffers
  //!  bc.serialize(A);  // forwards buffer into object A
  //!  bc.serialize(B);  // forwards buffer into object B
  //!  bc.serialize(C);  // forwards buffer into object C
  //!  bc.reset();       // Destroy buffers, bc ready for new transaction
  //! \endcode
  //! It is essential to load and forward data in the same order!! 
  //! A second interface uses the mpi::operator<<( BroadCast &, some type) to
  //! make this more simple. The above becomes:
  //! \code
  //!  mpi::BroadCast bc(mpi::main); // initializes a BroadCast objet named bc
  //!  bc << A << B << C  // sizes up object
  //!     << mpi::BroadCast::allocate // allocates buffers
  //!     << A << B << C  // loads buffers
  //!     << mpi::BroadCast::broadcast // broadcasts buffers
  //!     << A << B << C  // forwards buffers to A, B, C
  //!     << mpi::BroadCast::clear; // clears all and resets bc
  //! \endcode
  //! Whithin this code, only a few types can serialized, types::int,
  //! types::real, bool, and std::string. You can (and should) "overload" the
  //! template member function BroadCast::serialize() to include your own types
  //! and classes.
  //! \warning It is important that, when using BroadCast, containers on the
  //!          receiving-end processors should be empty (eg the container
  //!          should be empty on all nodes but the node which is
  //!          broadcasting). Otherwise, a buffer overflow error may occur when
  //!          copying a container to the buffer.
  //! \todo Make the problem of emptying containers on non-broadcasting node
  //!       transparent to the user somehow.
  class BroadCast : public CommBase 
  {
    //! \brief interface to this class! see class description for use
    //! \details handles both object serialization and broadcasting
    //! operation in one clear call.
    //! \param _this a BroadCast Object
    //! \param _type Make sure BroadCast::serialize() has been declared for your T_TYPE.
    template< class T_TYPE > friend BroadCast& operator<< ( BroadCast& _this, T_TYPE &_type );
    public:
      //! Defines all stages for BroadCast transactions
      enum t_stages { GETTING_SIZE,       //!< Size-up stage
                      COPYING_TO_HERE,    //!< Loading up stage
                      COPYING_FROM_HERE   //!< Forwarding stage
                    };
    protected:
      //! Defines all operations for BroadCast transactions
      enum t_operation { SIZEUP,     //!< Size-up 
                         ALLOCATE,   //!< Buffer allocation
                         TOHERE,     //!< Loading up buffers
                         BROADCAST,  //!< broadcasting buffers
                         FROMHERE,   //!< forwarding buffers
                         CLEAR,      //!< destory buffers, reinitialize object
                         NEXTSTAGE   //!< move to next stage
                       };
    public:
      //! allow calls to  BroadCast::CLEAR operation
      const static t_operation clear;    
      //! allow calls to BroadCast::NEXTSTAGE operation
      const static t_operation nextstage;
      //! allow calls to BroadCast::SIZEUP operation
      const static t_operation sizeup;
      //! allow calls to BroadCast::ALLOCATE operation
      const static t_operation allocate;
      //! allow calls to BroadCast::TOHERE operation
      const static t_operation tohere;
      //! allow calls to BroadCast::BROADCAST operation
      const static t_operation broadcast;
      //! allow calls to BroadCast::FROMHERE operation
      const static t_operation fromhere;

    protected:
      t_stages stage;  //!< Tracks which stage we are at
      bool keep_going; //!< upon error, shortcuts operator<< stuff

    public:
      //! Constructor, should not be used
      BroadCast() : CommBase(), stage(GETTING_SIZE), keep_going(true) {}
      //! Copy Constructor
      BroadCast( const BroadCast &_b ) : CommBase(_b), stage(_b.stage), keep_going(true) {}
      //! Initialization Constructor, should be used in conjunction with mpi::main
      BroadCast( const Base &_b ) : CommBase(_b), stage(GETTING_SIZE), keep_going(true) {}
      //! Destructor, destroy CommBase ineherited buffers
      virtual ~BroadCast() { destroy_buffers(); };

      //! \brief Allocate buffers
      //! \details allocate_buffers(types::t_unsigned _root) does just that. May
      //! become non-virtual in the future.
      //! \param _root optionally, specifies root process of data transaction
      bool allocate_buffers( types::t_unsigned _root = ROOT_NODE );

      //! \brief Allows to pass data to this class
      //! \details serialize() is used for sizing-up data, loading-up data, and then forwarding data
      //! It is defined only for a few basic types: you need to "overload" it for your own use
      //! \param _object object to be passed, make sure template function has been explicitely 
      //!        defined for T_OBJET!! Otherwise, you'll get a linkage error ;)
      template< class T_OBJECT > bool serialize( T_OBJECT& _object );
      //! \brief Serialize function for iterator ranges.
      //! \details Range must be valid from _first to _last, excluding _last.
      //! Furthermore a BroadCast::serialize() must be explicitely declared which
      //! takes T_ITERATOR::value_type as an argument.
      //! \param _first object to be serialized
      //! \param _last object in range, will \em not be  serialized
      template< class T_ITERATOR > bool serialize( T_ITERATOR _first, T_ITERATOR _last );
      //! Use this to go to next stage
      bool operator()( types::t_unsigned _root = ROOT_NODE );
      //! destorys buffers, and resets object to stage 0
      void reset() { destroy_buffers(); stage = GETTING_SIZE; keep_going = true; }

      //! returns which stage we are at, can be usefull in definition of serialize functions
      t_stages get_stage() const { return stage; } 

      //! \brief Serialize function for containers
      //! \details serialize_container() does just that. Keeps track of size of
      //! container. Upon loading, container will be  resized
      //! \param _cont to be serialized. Input and Output
      template< class T_CONTAINER >
      bool serialize_container ( T_CONTAINER &_cont );


      //! Point to Point transfer: Sends buffers to target. 
      void send_ptp( types::t_unsigned _target );
      //! Point to Point transfer: Receive buffers to target. 
      void receive_ptp( types::t_unsigned _source );

    protected:
      //! \brief Don't touch unless you know what you are doing
      //! \details allows operator<<(BroadCast &, T_TYPE&) to call both operations and to
      //! serialize objects
      template<class T_TYPE> bool operator_( T_TYPE &_type );
  };


  //! \brief Can be used to send gather data from all processes and then 
  //! distribute it to all processes.
  //! \details Can be used to send compound data from all processes in root
  //! process, and then send complete data back to all processes. It
  //! adopts the same six fold scheme as mpi::BroadCast :
  //!  - all variables should be sized up using BroadCast::serialize(), in
  //!      order to determine the size of the buffers. Any number of calls to 
  //!      BroadCast::serialize() can be made, with any type of argument for
  //!      which it has been explicitely defined,
  //!  - then the buffers should be allocated using AllGather::operator()(),
  //!  - a second series of calls to BroadCast::serialize should be made in
  //!      order to load up the buffers. 
  //!  - the data in the buffers should be gathered and distributed using
  //!      AllGather::operator()(),
  //!  - the data should be forwarded from the buffers to the serizalized
  //!      object using BroadCast::serialize(), in exactly the same order as
  //!      they were loaded in.
  //!  - a last optional call to AllGather::operator()() readies the object for
  //!      another data transaction.
  //!  .
  //! For instance, we have
  //! \code
  //!  mpi::AllGather ag(mpi::main); // initializes an AllGather objet named ag
  //!  ag.serialize(A);  // sizes up object A
  //!  ag.serialize(B);  // sizes up object B
  //!  ag.serialize(C);  // sizes up object C
  //!  ag.allocate_buffers(); // allocates buffers
  //!  ag.serialize(A);  // loads object A into buffers
  //!  ag.serialize(B);  // loads object B into buffers
  //!  ag.serialize(C);  // loads object C into buffers
  //!  ag(); // broadcasts buffers
  //!  ag.serialize(A);  // forwards buffer into object A
  //!  ag.serialize(B);  // forwards buffer into object B
  //!  ag.serialize(C);  // forwards buffer into object C
  //!  ag.resrt();       // Destroy buffers, ag ready for new transaction
  //! \endcode
  //! It is essential to load and forward data in the same order!! 
  //! No << interface has been designed at this point!
  //! Whithin this code, only a few types can serialized, types::int,
  //! types::real, bool, and std::string. You can (and should) "overload" the
  //! template member function BroadCast::serialize() to include your own types
  //! and classes.
  class AllGather : public BroadCast
  {
    //! \brief interface to this class! see class description for use
    //! \details handles both object serialization and broadcasting
    //! operation in one clear call.
    //! \param _this a BroadCast Object
    //! \param _type Make sure BroadCast::serialize() has been declared for your T_TYPE.
    template< class T_TYPE > friend AllGather& operator<< ( AllGather& _this, T_TYPE &_type );
    protected: 
      types::t_int *all_sizes; //!< receives gathered buffer sizes
    public:
      //! Constructor, don't use
      AllGather() : BroadCast(), all_sizes(NULL) {}
      //! Copy Constructor
      AllGather( const AllGather &_b ) : BroadCast(_b), all_sizes(NULL) {}
      //! Initialization Constructor, use in conljunction with mpi::main
      AllGather( const Base &_b ) : BroadCast(_b), all_sizes(NULL) {}
      //! Destructor, destroy CommBase buffers and size buffer
      virtual ~AllGather();
      //! \brief Allocate buffers
      //! \details allocate_buffers(types::t_unsigned _root) does just that. May
      //! become non-virtual in the future.
      //! \param _root optionally, specifies root process of data transaction
      bool allocate_buffers( types::t_unsigned _root = ROOT_NODE );
      //! Use this to go to next stage
      //! \param _root optionally, specifies root process of data transaction
      bool operator()( types::t_unsigned _root = ROOT_NODE );

    protected:
      //! \brief Don't touch unless you know what you are doing
      //! \details allows operator<<(BroadCast &, T_TYPE&) to call both operations and to
      //! serialize objects
      template<class T_TYPE> bool operator_( T_TYPE &_type );
  };

}

#include "comm.impl.h"
#endif
