//
//  Version: $Id$
//
#ifndef  _GRAPH_TOPOLOGY_H_
#define  _GRAPH_TOPOLOGY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef _MPI

#include <boost/mpi/graph_communicator.hpp>

#include <tinyxml/tinyxml.h>

#include <opt/types.h>
#include <opt/debug.h>
#include <print/stdout.h>
#include <mpi/mpi_object.h>

namespace GA
{
  class Topology; 

  namespace mpi 
  {
    //! \brief Graph topology.
    //! \details The topology and subsequent groups can be explained either
    //!          poetically or mathematically. One could imagine a farm headed
    //!          by a lonely single man. He believes himself intrusted by God
    //!          with the propagation of the \e Bos \e Taurus. He was granted a
    //!          number of bulls. And each bull a number of cows with which to
    //!          do the will of the Almighty. In his wisdom, the Master
    //!          Universe may have afflicted good-for-nothing farmhands to the
    //!          good man. More specifically, if there are \f$ n = p*h + 1\f$
    //!          processors. Processor +1 is the head boss directing \e p pools
    //!          of \e h processor. The processors of a pool work together to
    //!          evaluate a configuration. In any case, the topology consists
    //!          of an origin (the farmer) and \p rings. The origin is
    //!          connected to a single member of each ring. That's the bull.
    //!          The cows make up the other members of the ring. Rings
    //!          consisting of asingle bull are allowed (ring becomes a point),
    //!          as well as rings consisting of a cow and a bull ( ring becomes
    //!          a segment). 
    //
    //!          There are three "inputs" to the graph: \e n + 1 the number of
    //!          processes, \p the requested number of pools, \a _condition a
    //!          condition which \e h must fulfill (for instance, \t escan
    //!          requires that the number of waves is divisible by the number of
    //!          procs in a pool). Hence the farmhands. The algorithm strives
    //!          to find the best combination such that \b at \b most \e p pools
    //!          are created, and that these pools have a balanced number of
    //!          processor. It applies the condition (\a _condition in
    //!          Graph::Topology::init()) \e sine \e qua \e non.
    namespace Graph
    {
      struct t_Type
      {
        enum Process 
        {
          FARMER, //!< The head boss.
          FARMHAND, //!< Will do nothing.
          BULL,   //!< Takes orders from the farmer.
          COW     //!< Takes orders from one specific bull
        };
      };

      //! One possible condition for Topology::init()
      struct AlwaysTrue
      {
        //! Always returns true
        bool operator()( types::t_unsigned i ) { return true; }
      };

      //! Creates and contains the graph topology itself.
      class Topology 
      {
        friend class GA::Topology;
        protected: 
          //! The communicator for which processes will be graphed.
          boost::mpi::communicator *comm;
          //! If set, group of bull+cows of this process.
          boost::mpi::communicator pool_comm;
          //! If set, group of farmer+bulls of this process.
          boost::mpi::communicator head_comm;
          //! Graph communicator created by a call to MPI graph routine.
          boost::mpi::graph_communicator  graph_comm;
          //! The number of pools (e.g. herds) in the graph topology.
          types::t_unsigned pools;
          //! Type of the process
          t_Type::Process type;
      
        public:  
          //! \brief Constructor and Initializer.
          //! \details The MPI::WORLD_COMM is duplicated and the duplicated is
          //!          pointed to by ::mpi::Base::comm.
          Topology() : comm(::mpi::main), pool_comm(MPI::COMM_NULL), 
                       head_comm(MPI::COMM_NULL), graph_comm(MPI::COMM_NULL), 
                       pools(0), type(t_Type::FARMHAND)
             { comm = &MPI::COMM_WORLD.Clone(); }
          //! Copy Constructor.
          Topology   ( boost::mpi::communicator &_comm )
                   : comm(&_comm), pool_comm(MPI::COMM_NULL), 
                     head_comm(MPI::COMM_NULL), graph_comm(MPI::COMM_NULL), 
                     pools(0), type(t_Type::FARMHAND) {}
          //! Destructor
          ~Topology();
      
          //! \brief Creates the GA mpi topology and groups.
          template< class T_CONDITION > bool init( T_CONDITION _condition );
          //! Creates the mpi topology with the condition AlwaysTrue
          bool init () { return init( AlwaysTrue() ); }

      
          //! \brief sends random seeds to processes
          //! \details The seeds are resized to the number of pools+1. Any seed
          //!          equal to zero implies that seeding is done from the
          //!          current time.
          void reseed( std::vector< types::t_int > &_seeds );

          //! Loads the number of pools
          bool Load( const TiXmlElement &_node );

          //! Returns pointer to farmer+bull intracomm.
          mpi::boost::communicator* farmer_comm() { return &head_comm; }
          //! Returns pointer to bull/cows intracomm.
          mpi::boost::communicator* herd_comm() { return &pool_comm; }
          //! Sets mpi communicator and suffix.
          template< class T_EVALUATOR > void set_mpi( T_EVALUATOR &_eval );
      };

    } // namespace Graph


  } // namespace mpi
} // namespace GA

#include "graph.impl.h"

#endif
#endif
