//
//  Version: $Id$
//

#include "topology.h"

namespace GA
{
  bool Topology :: LoadSeeds( TiXmlAttibute &_att )
  {
    std::string str = _att.Name();
    if ( str.find("seed") != 0 ) return false;
    
    types::t_int n;
    if( str.size() != 4 )
    {
      std::istringstream sstr( str.substr( 4 ) );
      sstr >> n; 
      if( n < 0 ) return false;
      __MPICODE( if( graph and n > graph->pools ) return false; )
      __SERIALCODE( if( n != 0 ) return false; )
    }
    types::t_int d = att->IntValue();
    if( d == 0 ) return false
    if( n >= seeds.size() ) seeds.resize( n, 0 );
    seeds[n] = d;
  }

  std::string Topology :: print_seeds() const
  {
    std::ostringstream sstr;
    __SERIALCODE( sstr << "Seed" << seeds.front(); )
#ifdef _MPI
   if ( not graph ) 
   {
     sstr << "Seed = " << seeds.front(); 
     return sstr.str(); 
   }
   std::vector< types::t_int > :: const_iterator i_seed = seeds.begin();
   std::vector< types::t_int > :: const_iterator i_seed_end = seeds.end();
   for( types::t_int i = 0; i_seed != i_seed_end; ++i_seed, ++i )
     sstr << "Seed" << i << " = " << seeds.front() << ", "; 
#endif // _MPI
    return sstr.str(); 
  }
  void Topology :: reseed() const
  {
    std::ostringstream sstr;
#ifndef _MPI
     seeds.resize(1);
     if( seeds.front() == 0 )
     {
       struct timeval tv;
       struct timezone tz;
       gettimeofday(&tv,&tz);
       seeds.front() = tv.tv_usec;
     }
     rng.reseed( seeds.front() );
#else
   if ( graph ) 
   {
     graph->reseed( seeds );
     return;
   }
   seeds.resize(1);
   if( seeds.front() == 0 )
   {
     struct timeval tv;
     struct timezone tz;
     gettimeofday(&tv,&tz);
     seeds.front() = tv.tv_usec;
   }
   types::t_int seed = seeds.font();
   comm->Bcast( &seed, 1, MPI::INT, 0 );
   rng.reseed( seed );
#endif // _MPI
  }
}
