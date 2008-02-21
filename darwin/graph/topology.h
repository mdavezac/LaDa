//
//  Version: $Id$
//
#ifndef  _GRAPH_TOPOLOGY_H_
#define  _GRAPH_TOPOLOGY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef _MPI

#include <opt/types.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>

namespace GA
{
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
      //! One possible condition for Topology::init()
      struct AlwaysTrue
      {
        //! Always returns true
        bool operator()( types::t_unsigned i ) { return true; }
      };

      //! Creates and contains the graph topology itself.
      class Topology : private ::mpi::Base
      {
        protected: 
          //! If set, group of bull+cows of this process.
          MPI::IntraComm *herd_comm;
          //! If set, group of farmer+bulls of this process.
          MPI::IntraComm *farmer_comm;
          //! Graph communicator created by a call to MPI graph routine.
          MPI::GraphComm  *graph_comm;
          //! The number of pools (e.g. herds) in the graph topology.
          types::t_unsigned pools;
          //! Type of the process
          enum t_Type 
          {
            FARMER, //!< The head boss.
            FARMHANDS, //!< Will do nothing.
            BULL,   //!< Takes orders from the farmer.
            COW     //!< Takes orders from one specific bull
          }
          t_Type type;
      
        public:  
          //! \brief Constructor and Initializer.
          //! \details The MPI::WORLD_COMM is duplicated and the duplicated is
          //!          pointed to by ::mpi::Base::comm.
          Topology() : ::mpi::Base::Base(), herd_comm(NULL), farmer_comm(NULL),
                       graph_comm(NULL), pools(0), type(FARMHANDS)
             { comm = &_comm->Clone(); }
          //! Copy Constructor.
          Topology   ( MPI::Comm &_comm )
                   : ::mpi::Base::Base(_comm), herd_comm(NULL), farmer_comm(NULL),
                     graph_comm(NULL), pools(0), type(FARMHANDS)
            { comm = &_comm->Clone(); }
          //! Destructor
          ~Topology();
      
          //! \brief Creates the GA mpi topology and groups.
          template< class T_CONDITION > void init( T_CONDITION &_condition );
          //! Creates the mpi topology with the condition AlwaysTrue
          void init () { init( AlwaysTrue() ); }

      
          //! \brief sends random seeds to processes
          //! \details The seeds are resized to the number of pools+1. Any seed
          //!          equal to zero implies that seeding is done from the
          //!          current time.
          void reseed( std::vector< types::t_int > &_seeds );

          //! Loads the number of pools
          bool Load( const TiXmlElement &_node );
      
      };

    } // namespace Graph


  } // namespace mpi
} // namespace GA

#include "topology.impl.h"

#endif
#endif
