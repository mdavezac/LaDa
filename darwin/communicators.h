//
//  Version: $Id$
//
#ifndef  _DARWIN_COMMUNICATORS_H_
#define  _DARWIN_COMMUNICATORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>

namespace GA
{
  namespace mpi 
  {

    class Graph : ::mpi::Base
    {
      protected: 
        MPI::Comm *group_comm;
        MPI::GraphComm  *graph_comm;
        types::t_unsigned pools;
        t_Condition& condition;
        //! Type of the process
        enum t_Type 
        {
          COLLECTIVIST, //!<  Equivalent to pool == 1.
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
        Graph() : ::mpi::Base::Base(), group_comm(NULL),
                  graph_comm(NULL), pools(0), type(FARMHANDS) { comm = &_comm->Clone(); }
        Graph   ( MPI::Comm &_comm )
              : ::mpi::Base::Base(_comm), group_comm(NULL),
                graph_comm(NULL), pools(0), type(FARMHANDS) { comm = &_comm->Clone(); }
        ~Graph();
        
        MPI::Comm& Comm() { return *group_comm; }
    };

    template< class T_CONDITION >
    Graph* graph_init( MPI::Comm &_comm , t_Condition &)


  } // namespace mpi
} // namespace GA

#endif
