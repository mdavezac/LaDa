//
//  Version: $Id$
//

#include <sys/time.h>
#include <eo/utils/eoRNG.h>

#include "graph.h"

namespace GA
{
  namespace mpi 
  {
    namespace Graph
    {
      Topology :: ~Topology()
      {
        if( graph_comm != *comm and graph_comm != MPI::COMM_NULL )
          graph_comm.Free();
        if( comm ) comm->Free();
        if( head_comm != MPI::COMM_NULL ) head_comm.Free();
        if( pool_comm != MPI::COMM_NULL ) pool_comm.Free();
        graph_comm = NULL;
        comm = NULL;
        pool_comm = NULL;
        head_comm = NULL;
      }
     
      void Topology :: reseed( std::vector< types::t_int > &_seeds )
      {
        if( type == t_Type::FARMER )
        {
          _seeds.resize( pools + 1 );
          std::vector< types::t_int > :: iterator i_seed = _seeds.begin();
          std::vector< types::t_int > :: iterator i_seed_end = _seeds.end();
          for(; i_seed != i_seed_end; ++i_seed )
            if( *i_seed == 0 )
            {
              struct timeval tv;
              struct timezone tz;
              gettimeofday(&tv,&tz);
              *i_seed = tv.tv_usec;
            }
          i_seed = _seeds.begin();
          rng.reseed( *i_seed ); ++i_seed;
          for( types::t_unsigned i = 0; i_seed != i_seed_end; ++i_seed, ++i)
            head_comm.Send( &(*i_seed), 1, MPI::INT, i, 666 );
        }
        else if( type == t_Type::BULL )
        {
          types::t_int seed;
          comm->Recv( &seed, 1, MPI::INT, 0, 666 );
          pool_comm.Bcast( &seed, 1, MPI::INT, 0 );
          rng.reseed( seed );
        }
        else if ( type == t_Type::COW )
        {
          types::t_int seed;
          comm->Bcast( &seed, 1, MPI::INT, 0 );
          rng.reseed( seed );
        }
      }
     
      bool Topology :: Load ( const TiXmlElement &_node )
      {
        const TiXmlElement *child = _node.FirstChildElement( "MPI" );
        if ( not ( child and child->Attribute("pools") ) ) return false;
        types::t_int i;
        child->Attribute( "pools", &i );
        if ( i < 1 ) return false;
        pools = std::abs(i);
        return true;
      }

    } // namespace Graph


  } // namespace mpi
} // namespace GA



