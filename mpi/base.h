//
//  Version: $Id$
//
#ifdef _MPI
#ifndef _MPI_BASE_H_
#define _MPI_BASE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>

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

  //! Simple %Base class to all other helper classes.
  class Base
  {
    protected:
      types::t_int this_rank; //!< rank of current process.
      types::t_int nproc; //!< size of mpi pool
      MPI::Comm *comm;

    public:
       //! Constructor and Initializer
      Base () : this_rank(0), nproc(1), comm(&MPI::COMM_WORLD)  {}
       //! Constructor and Initializer
      Base   ( MPI::Comm &_comm ) 
           : this_rank( _comm.rank() ), nproc( _comm.size() ),
             comm(&_comm)  {}
      //! Copy Constructor
      Base   ( const Base &_c)  
           : this_rank( _c.this_rank), nproc( _c.nproc), comm(_c.comm) {}
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
      void barrier() const  { comm->Barrier(); }

      //! \brief sums an unsigned integer across all processes
      //
      //! \param _in  unsigned integer argument to be summed across all procs,
      //!             is input and output
      types::t_unsigned all_sum_all( types::t_unsigned &_in) const;
      //! \brief sums an integer across all processes
      //
      //! \param _in  integer argument to be summed across all procs,
      //!             is input and output
      types::t_int all_sum_all( types::t_int &_in) const;
      //! \brief sums a real value across all processes
      //
      //! \param _in  real value to be summed across all procs,
      //!             is input and output
      types::t_real all_sum_all( types::t_real &_in) const;
      //! \brief sums bool value across all processes, eg returns true if
      //argument is true across all processes
      //
      //! \param _bool  bool value to be tested for true across all procs,
      //!             is input and output
      bool all_sum_all( bool &_bool) const;
      //! \brief returns true if argument is true on at least one process.
      //
      //! \param _bool  bool value to be tested for true across all procs,
      //!             is input and output
      bool all_or_all( bool &_bool) const;
      //! \brief returns true if argument is either true across all processes
      //         or none of the processes
      //
      //! \param _bool  bool value to be tested for true across all procs,
      //!             is input and output
      bool all_xor_all( bool &_bool) const;
      //! Sets the commiunication handle
      void set_comm( Mpi::Comm *_comm )
       { comm = _comm; this_rank = comm->rank(); nproc = comm->size(); }
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
        this_rank = comm->Get_rank();
        nproc = comm->Get_size();
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

  inline Base :: set_comm( Mpi::Comm *_comm )
  { 
    __ASSERT( _comm, "Communicator is null,\n" )
    comm = _comm;
    this_rank = comm->rank(); 
    nproc = comm->size(); 
  }
  inline types::t_unsigned Base::all_sum_all( types::t_unsigned &_in) const
  {
    types::t_unsigned out;
    comm->Allreduce( &_in, &out, 1, UNSIGNED, MPI::SUM ); 
    _in = out;
    return out;
  }

  inline types::t_int Base::all_sum_all( types::t_int &_in) const
  {
    types::t_int out;
    comm->Allreduce( &_in, &out, 1, INT, MPI::SUM ); 
    _in = out;
    return out;
  }

  inline types::t_real Base::all_sum_all( types::t_real &_in) const
  {
    types::t_real out;
    comm->Allreduce( &_in, &out, 1, REAL, MPI::SUM ); 
    _in = out;
    return out;
  }

  inline bool Base::all_sum_all( bool &_bool) const
  {
    types::t_int out, in = _bool ? 1 : 0;
    comm->Allreduce( &in, &out, 1, UNSIGNED, MPI::SUM ); 
    _bool = ( out == nproc );
    return _bool;
  }

  inline bool Base::all_or_all( bool &_bool) const
  {
    types::t_int out, in = _bool ? 1 : 0;
    comm->Allreduce( &in, &out, 1, UNSIGNED, MPI::SUM ); 
    _bool = ( out != 0 );
    return _bool;
  }

  inline bool Base::all_xor_all( bool &_bool) const
  {
    types::t_int out, in = _bool ? 1 : 0;
    comm->Allreduce( &in, &out, 1, UNSIGNED, MPI::SUM ); 
    _bool = ( out == 0 or out == nproc );
    return _bool;
  }

  inline void InitDestroy::operator()( int _argc, char **_argv )
  { 
    if ( MPI::Is_initialized() or finalized )
      return;
    MPI::Init( _argc, _argv );
    this_rank = comm->Get_rank();
    nproc = comm->Get_size();
  }
  
  inline InitDestroy::~InitDestroy()
  { 
    if ( ( not MPI::Is_initialized() ) or finalized )
      return;
    MPI::Finalize();
    finalized = true;
  }

} // namespace mpi
/*@}*/
#endif
