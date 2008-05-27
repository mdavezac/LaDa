//
//  Version: $Id$
//

#include <sys/time.h>
#include <eo/utils/eoRNG.h>

#include "graph.h"

#ifdef _MPI
#include <boost/mpi/collectives.hpp>
#endif

namespace GA
{
  namespace mpi 
  {
    namespace Graph
    {
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
          for( types::t_unsigned i = 1; i_seed != i_seed_end; ++i_seed, ++i)
            head_comm.send( i, 666, *i_seed);
        }
        else if( type == t_Type::BULL )
        {
          types::t_int seed;
          head_comm.recv( 0, 666, seed);
          boost::mpi::broadcast( pool_comm, seed, 0 );
          rng.reseed( seed );
        }
        else if ( type == t_Type::COW )
        {
          types::t_int seed;
          boost::mpi::broadcast( pool_comm, seed, 0 );
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



