//
//  Version: $Id$
//
#ifndef _MPI_OBJECT_H_
#define _MPI_OBJECT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdexcept>       // std::runtime_error
#include <iostream>
#include <math.h>
#ifdef _MPI
#include <openmpi/mpi.h>
//<mpi2c++/mpi++.h>
#endif

#include "opt/types.h"

/** \ingroup MPI
 * @{*/
//!\brief mpi namespace should contain all mpi related helper classes. 
/*!       Currently works with openmpi implementation only. 
  
     You will need mpich header files and libraries
     At present the following helper classes are defined:
      - mpi::Base : base class to all others, it contains functions to return rank and
              size of the current process, as well as simple all_sum_all
              functions for single variables.
      - mpi::InitDestroy : helper class for creating and destroying mpi. should
                     probably not be derived from.
      - mpi::CommBase : Base for data communication helper classes. You can use this class
                  to further create comm helpers.
         - mpi::BroadCast : helper for broadcasting, eg sending data from one process
                            to all processes
         - mpi::AllGather : helper for gathering and distributing data, eg
                            sending data from all processes to one processes and
                            back to all processes. In the end all process should
                            have access to the same data.
*/
namespace mpi
{
#define REAL MPI::DOUBLE //!< renames the openmpi macro for double values
#define INT MPI::INT     //!< renames the openmpi macro for integer values
#define CHAR MPI::CHAR   //!< renames the openmpi macro for character values
//! renames the openmpi macro for unsigned integer values
#define UNSIGNED MPI::UNSIGNED

//! one root to rule them all, a priori
extern const types::t_int ROOT_NODE;

#ifdef _MPI
  //! Simple %Base class to all other helper classes.
  class Base
  {
    protected:
      types::t_int this_rank; //!< rank of current process.
      types::t_int nproc; //!< size of mpi pool

    public:
      Base () : this_rank(0), nproc(1) {} //<! constructor 
      //! Copy Constructor
      Base   ( const Base &_c)  
           : this_rank( _c.this_rank), nproc( _c.nproc) {}
      virtual ~Base() {} //!< destructor

      //! returns rank of current process
      types::t_int rank() const { return this_rank; } 
      //! returns size of mpi pool
      types::t_int size() const { return nproc; }
      
      //! returns true if node is root, as defined by mpi::ROOT_NODE
      bool is_root_node() const { return this_rank == ROOT_NODE; }

      //!\brief simple barrier function.
      //
      //! Processes stop when mpi::barrier() const is called until all process
      //! have reached mpi::barrier() const
      void barrier() const
        { MPI::COMM_WORLD.Barrier(); }

      //! \brief sums an unsigned integer across all processes
      //
      //! \param _in  unsigned integer argument to be summed across all procs,
      //!             is input and output
      types::t_unsigned all_sum_all( types::t_unsigned &_in) const
      {
        types::t_unsigned out;
        MPI::COMM_WORLD.Allreduce( &_in, &out, 1, UNSIGNED, MPI::SUM ); 
        _in = out;
        return out;
      }

      //! \brief sums an integer across all processes
      //
      //! \param _in  integer argument to be summed across all procs,
      //!             is input and output
      types::t_int all_sum_all( types::t_int &_in) const
      {
        types::t_int out;
        MPI::COMM_WORLD.Allreduce( &_in, &out, 1, INT, MPI::SUM ); 
        _in = out;
        return out;
      }

      //! \brief sums a real value across all processes
      //
      //! \param _in  real value to be summed across all procs,
      //!             is input and output
      types::t_real all_sum_all( types::t_real &_in) const
      {
        types::t_real out;
        MPI::COMM_WORLD.Allreduce( &_in, &out, 1, REAL, MPI::SUM ); 
        _in = out;
        return out;
      }

      //! \brief sums bool value across all processes, eg returns true if
      //argument is true across all processes
      //
      //! \param _bool  bool value to be tested for true across all procs,
      //!             is input and output
      bool all_sum_all( bool &_bool) const
      {
        types::t_int out, in = _bool ? 1 : 0;
        MPI::COMM_WORLD.Allreduce( &in, &out, 1, UNSIGNED, MPI::SUM ); 
        _bool = ( out == nproc );
        return _bool;
      }

      //! \brief returns true if argument is true on at least one process.
      //
      //! \param _bool  bool value to be tested for true across all procs,
      //!             is input and output
      bool all_or_all( bool &_bool) const
      {
        types::t_int out, in = _bool ? 1 : 0;
        MPI::COMM_WORLD.Allreduce( &in, &out, 1, UNSIGNED, MPI::SUM ); 
        _bool = ( out != 0 );
        return _bool;
      }

      //! \brief returns true if argument is either true across all processes
      //         or none of the processes
      //
      //! \param _bool  bool value to be tested for true across all procs,
      //!             is input and output
      bool all_xor_all( bool &_bool) const
      {
        types::t_int out, in = _bool ? 1 : 0;
        MPI::COMM_WORLD.Allreduce( &in, &out, 1, UNSIGNED, MPI::SUM ); 
        _bool = ( out == 0 or out == nproc );
        return _bool;
      }
  };

  //! \brief Handles creation and destruction of mpi aspect
  //! \details A global mpi::InitDestroy object is created, mpi::main, which should be
  //! called once, preferably at the start of the main() function. Simply
  //! place call mpi::InitDestroy::operator()(int, char **) at the start of any
  //! main(int, char**) with the corresponding arguments and mpi will be
  //! initialized. mpi are finalized upon the destruction of mpi::main
  class InitDestroy : public Base
  {
    protected:
      bool finalized; //!< checks wether mpi has already been finalized
    public:
      //! constructor, initailizes mpi::InitDestroy::finalized to false
      InitDestroy () : Base (), finalized(false) {} 
      //! \brief should be called once to enable mpi 
      //! \details Enables mpi and initializes Base::this_rank and Base::nproc
      //!    \param _argc see main(int, char**)
      //!    \param _argv see main(int, char**)
      void operator()( int _argc, char **_argv )
      { 
        if ( MPI::Is_initialized() or finalized )
          return;
        MPI::Init( _argc, _argv );
        this_rank = MPI::COMM_WORLD.Get_rank();
        nproc = MPI::COMM_WORLD.Get_size();
      }
      //! \brief Destructor, disables mpi calls
      //! \details Disables mpi calls. After this function is called, no other mpi calls
      //! should be made.
      virtual ~InitDestroy()
      { 
        if ( ( not MPI::Is_initialized() ) or finalized )
          return;
        MPI::Finalize();
        finalized = true;
      }
  };

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
      CommBase () : Base(), int_buff(NULL), end_int_buff(NULL), cur_int_buff(NULL),
                    char_buff(NULL), end_char_buff(NULL), cur_char_buff(NULL),
                    real_buff(NULL), end_real_buff(NULL), cur_real_buff(NULL)
      { 
        buffer_size[0] = 0;
        buffer_size[1] = 0;
        buffer_size[2] = 0;
      }
      //! \brief initialization Constructor, copies only mpi::Base members.
      //! \details Should be used in conjunction with mpi::main to initialize any data
      //! communication helper class. Only mpi::Base members are copied,
      //! Members from CommBase are set to 0 or NULL. 
      CommBase   ( const Base &_c ) 
               : Base(_c), int_buff(NULL), end_int_buff(NULL), cur_int_buff(NULL),
                 char_buff(NULL), end_char_buff(NULL), cur_char_buff(NULL),
                 real_buff(NULL), end_real_buff(NULL), cur_real_buff(NULL)
      { 
        buffer_size[0] = 0;
        buffer_size[1] = 0;
        buffer_size[2] = 0;
      }
      //! \brief true Copy Constructor
      CommBase   ( const CommBase &_c) 
               : Base( _c ),
                 int_buff(_c.int_buff), end_int_buff(_c.end_int_buff), cur_int_buff(_c.cur_int_buff),
                 char_buff(_c.char_buff), end_char_buff(_c.end_char_buff), cur_char_buff(_c.cur_char_buff),
                 real_buff(_c.real_buff), end_real_buff(_c.end_real_buff), cur_real_buff(_c.cur_real_buff)
      {
        buffer_size[0] = _c.buffer_size[0];
        buffer_size[1] = _c.buffer_size[1];
        buffer_size[2] = _c.buffer_size[2];
      }
      //! \brief Destructor. Does nothing ;)
      //! \details Destructor leaves it to derived classes to destroy buffers!!
      virtual ~CommBase() {}

      //! Destroys all buffers if they exist, then assigns all members to 0 or NULL;
      void destroy_buffers()
      {
        if ( int_buff ) delete[] int_buff;
        if ( char_buff ) delete[] char_buff;
        if ( real_buff ) delete[] real_buff;
        int_buff = NULL;  end_int_buff = NULL;  cur_int_buff = NULL;  buffer_size[0] = 0;  
        real_buff = NULL; end_real_buff = NULL; cur_real_buff = NULL; buffer_size[1] = 0; 
        char_buff = NULL; end_char_buff = NULL; cur_char_buff = NULL; buffer_size[2] = 0; 
      }
     
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
      void reset()
        { destroy_buffers(); stage = GETTING_SIZE; keep_going = true; }

      //! returns which stage we are at, can be usefull in definition of serialize functions
      t_stages get_stage() const { return stage; } 

      //! \brief Serialize function for containers
      //! \details serialize_container() does just that. Keeps track of size of
      //! container. Upon loading, container will be  resized
      //! \param _cont to be serialized. Input and Output
      template< class T_CONTAINER >
      bool serialize_container ( T_CONTAINER &_cont )
      {
        types::t_int n = _cont.size();
        if ( not serialize( n ) ) return false;
        if( stage == COPYING_FROM_HERE )
          _cont.resize(n);
        typename T_CONTAINER :: iterator i_ob = _cont.begin();
        typename T_CONTAINER :: iterator i_ob_end = _cont.end();
        for(; i_ob != i_ob_end; ++i_ob )
          if( not serialize( *i_ob ) ) return false;
      
        return true;
      }
    protected:
      //! \brief Don't touch unless you know what you are doing
      //! \details allows operator<<(BroadCast &, T_TYPE&) to call both operations and to
      //! serialize objects
      template<class T_TYPE> bool operator_( T_TYPE &_type );
  };

  template<class T_TYPE> inline bool BroadCast :: operator_( T_TYPE &_type )
  {
    if ( not keep_going ) return true;
    try
    {
      if (not serialize( _type ) )
        throw std::runtime_error( "Error while BroadCasting \n" );
    }
    catch( std::exception &e )
    {
      std::cerr << "Caught error while running lada" << std::endl
                << e.what();
      keep_going = false;
      return false;
    }
    return true;
  }
  //! Takes care of operation calls for operator<<( BroadCast&, T_TYPE&)
  template<> inline bool BroadCast::operator_<const BroadCast::t_operation>( const t_operation &_op )
  {
    if ( not keep_going ) return true;
    if      ( _op == SIZEUP or _op == CLEAR )   reset();
    else if ( _op == ALLOCATE ) allocate_buffers(); 
    else if ( _op == BROADCAST ) operator()(); 
    else if ( _op == NEXTSTAGE )
    {
      switch( stage )
      {
        case GETTING_SIZE: allocate_buffers(); break;
        case COPYING_TO_HERE: operator()(); break;
        case COPYING_FROM_HERE: reset(); break;
        default: break;
      }
    }
    return true; 
  }

  template< class T_TYPE > inline BroadCast& operator<< ( BroadCast& _this, T_TYPE &_type )
    { _this.operator_( _type ); return _this; }

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
      virtual ~AllGather() 
      {
        if(all_sizes) delete[] all_sizes;
        all_sizes = NULL;
        destroy_buffers();
      };

      //! \brief Allocate buffers
      //! \details allocate_buffers(types::t_unsigned _root) does just that. May
      //! become non-virtual in the future.
      //! \param _root optionally, specifies root process of data transaction
      bool allocate_buffers( types::t_unsigned _root = ROOT_NODE );
      //! Use this to go to next stage
      //! \param _root optionally, specifies root process of data transaction
      bool operator()( types::t_unsigned _root = ROOT_NODE );
  };


  //! \brief Return an iterator so each processor gets to screen different bits of a container
  //! \details Each proc starts at a different point in the array, and ends right were the next begins
  //! \sa end(T_CONTAINER&), begin(T_CONTAINER&), begin(const T_CONTAINER&), end( const T_CONTAINER&)
  //! \param _cont a container which should be the same across all processors
  template<class T_CONTAINER>
  inline typename T_CONTAINER::iterator begin( T_CONTAINER &_cont )
  {
    extern InitDestroy main;
    types::t_real n =   (types::t_real) main.rank() 
                      * (types::t_real) _cont.size()
                      / (types::t_real) main.size();
    return _cont.begin() + ( types::t_unsigned ) floor(n);
  }
  //! \details Returns an iterator so each processor gets to screen different bits of a container
  //! \details Each proc starts at a different point in the array, and ends right were the next begins
  //! \sa end(T_CONTAINER&), begin(T_CONTAINER&), begin(const T_CONTAINER&), end( const T_CONTAINER&)
  //! \param _cont a container which should be the same across all processors
  template<class T_CONTAINER>
  inline typename T_CONTAINER::iterator end( T_CONTAINER &_cont )
  {
    extern InitDestroy main;
    types::t_real n =   (types::t_real) ( main.rank() + 1 )
                      * (types::t_real) _cont.size()
                      / (types::t_real) main.size();
    return _cont.begin() + ( types::t_unsigned ) floor(n);
  }
  //! \brief Returns an constant iterator so each processor gets to screen different
  //! bits of a container
  //! \details Each proc starts at a different point in the array, and ends right were the next begins
  //! \sa end(T_CONTAINER&), begin(T_CONTAINER&), begin(const T_CONTAINER&), end( const T_CONTAINER&)
  //! \param _cont a container which should be the same across all processors
  template<class T_CONTAINER>
  inline typename T_CONTAINER::const_iterator begin( const T_CONTAINER &_cont )
  {
    extern InitDestroy main;
    types::t_real n =   (types::t_real) main.rank() 
                      * (types::t_real) _cont.size()
                      / (types::t_real) main.size();
    return _cont.begin() + ( types::t_unsigned ) floor(n);
  }
  //! \brief Returns an constant iterator so each processor gets to screen different
  //! bits of a container
  //! \details Each proc starts at a different point in the array, and ends right were the next begins
  //! \sa end(T_CONTAINER&), begin(T_CONTAINER&), begin(const T_CONTAINER&), end( const T_CONTAINER&)
  //! \param _cont a container which should be the same across all processors
  template<class T_CONTAINER>
  inline typename T_CONTAINER::const_iterator end( const T_CONTAINER &_cont )
  {
    extern InitDestroy main;
    types::t_real n =   (types::t_real) ( main.rank() + 1 )
                      * (types::t_real) _cont.size()
                      / (types::t_real) main.size();
    return _cont.begin() + ( types::t_unsigned ) floor(n);
  }
#else

  class Base
  {
     public:
       Base () {};
       Base ( const Base &_c) {};
       virtual ~Base() {}
 
       types::t_int rank() const { return 0; }
       types::t_int size() const { return 1; }
       
       void barrier() const {};
       template< class T_TYPE >
       T_TYPE all_sum_all( T_TYPE &_un ) const 
         { return _un; };
       bool all_sum_all( bool &_bool ) const 
         { return _bool; };
       bool all_or_all( bool &_bool ) const 
         { return _bool; };
       bool all_xor_all( bool &_bool ) const 
         { return _bool; };
  };

  class CommBase : virtual public Base
  {
    public:
      CommBase () {};
      CommBase ( const CommBase &_c) {};
      CommBase ( const Base &_c) {};
      virtual ~CommBase() {}

      void destroy_buffers() {};
  };
 
  class InitDestroy : virtual public Base
  {
     public:
       InitDestroy () : Base() {};
       void operator()( int, char ** ) {};
       virtual ~InitDestroy() {};
  }; 
 
  class BroadCast : virtual public CommBase 
  {
    protected:
      enum t_stages { GETTING_SIZE, COPYING_TO_HERE, COPYING_FROM_HERE };
 
    public:
      BroadCast () : CommBase () {}
      BroadCast ( const BroadCast &_b ) {}
      BroadCast ( const Base &_b )  {}
      virtual ~BroadCast() {};
 
      
      bool allocate_buffers( types::t_unsigned _root = ROOT_NODE) { return true; };
 
      template< class T_OBJECT > bool serialize( T_OBJECT _object )
        { return true; }
      template< class T_ITERATOR > bool serialize( T_ITERATOR _first, T_ITERATOR _last )
        { return true; }
      template< class T_CONTAINER > bool serialize_container ( T_CONTAINER &_cont )
        { return true; }
      bool operator()( types::t_unsigned _root = ROOT_NODE ) { return true; };
      void reset() {};
  };
  class AllGather : public BroadCast
  {
    public:
      AllGather() : BroadCast() {}
      AllGather( const AllGather &_b ) : BroadCast(_b) {}
      AllGather( const Base &_b ) : BroadCast(_b) {}
      virtual ~AllGather() {}
  };
 
  template<class T_CONTAINER>
  inline typename T_CONTAINER::iterator begin( T_CONTAINER &_cont )
    { return _cont.begin(); } 
  template<class T_CONTAINER>
  inline typename T_CONTAINER::iterator end( T_CONTAINER &_cont )
    { return _cont.end(); }

#endif


 //! \brief Main object for creation and destruction of mpich, as well as for intialization
 //! of transaction helpers
 extern InitDestroy main;


    

}
/*@}*/

#endif
