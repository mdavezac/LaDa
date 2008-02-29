//
//  Version: $Id$
//

namespace GA
{
  namespace mpi 
  {
    namespace Graph
    {
      Topology :: ~Topology()
      {
        if( graph_comm and *comm != MPI::COMM_NULL and graph_comm != comm )
          graph_comm->Free();
        if( comm and *comm != MPI::COMM_NULL ) comm->Free();
        if( head_comm and head_comm != MPI::COMM_NULL ) head_comm->Free();
        if( pool_comm and pool_comm != MPI::COMM_NULL ) pool_comm->Free();
        graph_comm = NULL;
        comm = NULL;
        pool_comm = NULL;
        head_comm = NULL;
      }
     
      void Topology :: reseed( std::vector< types::t_int > &_seeds )
      {
        if( type == FARMER )
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
            head_comm->Send( &(*i_seed), 1, MPI::INT, i );
        }
        else if( type == BULL )
        {
          types::t_int seed;
          comm->Recv( &seed, 1, MPI::INT, 0 );
          cowcomm->Bcast( &seed, 1, MPI::INT, 0 );
          rng.reseed( seed );
        }
        else ( type == COW )
        {
          types::t_int seed;
          comm->Bcast( &seed, 1, MPI::INT, 0 );
          rng.reseed( seed );
        }
      }
     
      bool Topology :: Load ( const TiXmlElement &_node )
      {
        const TiXmlElement *child = _node.FirstChildElement( "MPI" );
        if ( not child->Attribute( "pools" ) ) return false;
        types::t_int i;
        child->Attribute( "pools", &i );
        if ( i < 1 ) return false;
        pools = std::abs(i);
        return true;
      }

    } // namespace Graph


  } // namespace mpi
} // namespace GA



