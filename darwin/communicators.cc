//
//  Version: $Id$
//

namespace GA
{
  namespace mpi 
  {
    Graph<T_CONDITION>::~Graph()
    { 
      if( graph_comm or comm or herd_comm or farmer_com )
        std::cerr << "Call to GA::mpi::Graph::destroy() is missing.\n";
    }
    Graph :: destroy()
    {
      if( graph_comm and *comm != MPI::COMM_NULL and graph_comm != comm )
        graph_comm->Free();
      if( comm and *comm != MPI::COMM_NULL ) comm->Free();
      if( farmer_comm and farmer_comm != MPI::COMM_NULL ) farmer_comm->Free();
      if( herd_comm and herd_comm != MPI::COMM_NULL ) herd_comm->Free();
      graph_comm = NULL;
      comm = NULL;
      herd_comm = NULL;
      farmer_comm = NULL;
    }


  }
}



