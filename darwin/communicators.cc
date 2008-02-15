//
//  Version: $Id$
//

namespace GA
{
  namespace mpi 
  {
    template< class T_CONDITION >
    Graph<T_CONDITION>::~Graph()
    { 
      if( comm ) comm->Free();
      if( graph_comm ) graph_comm->Free();
      if( group_comm ) group_comm->Free();
    }

    template< class T_CONDITION >
    void Graph<T_CONDITION>::init_collectivists()
    {
      group_comm = comm;
      types::t_unsigned per_pool = ::mpi::Base::size();
      while ( per_pool and ( not _condition( per_pool) ) ) --per_pool;
      if( per_pool == 0 )
        __THROW_ERROR(    "Could not satisfy the mandatory condition\n"
                       << "regarding the number of procs despite \n"
                       << "only using one pool\n" )
      types::t_unsigned nb = ::mpi::Base::size() - per_pool;
      Print::out << "MPI Process Organization: \n"
                 << "  _ Single pool of " << per_pool << " processes.\n"
                 << "  _ " << nb << " donothings.\n";
      type = COLLECTIVIST;
      if( not nb ) return;
      // distributes farmhands all over the spectrum;
      nb = ::mpi::Base::size() / nb;
      types::t_int color = (rank() + 1) % nb ? 0: MPI::UNDEFINED;
      type = color == MPI::UNDEFINED ? FARMHANDS: COLLECTIVIST;
      group_comm = &comm->Split( color, rank() );
      return;
    }

    template< class T_CONDITION >
    void Graph<T_CONDITION>::init()
    {
      //! One single pool
      if( pools <= 1 ) init_collectivists( _condition );
      types::t_unsigned per_pool = ( ::mpi::Base::size() - 1) / pools;
      types::t_unsigned leftovers = ( ::mpi::Base::size() - 1) % pools;

      while ( per_pool and ( not _condition( per_pool) ) ) --per_pool;

      //! Could not satisfy condition
      if( per_pool == 0 ) 
      {
        std::cerr << "Could not create " << pools << " pools of processors.\n"
                  << "Will try to create " << pools - 1 << " pools instead. \n";
        --pool;
        init( _condition );
        return;
      }

      typedef std::pair<types::t_unsigned, types::t_unsigned > t_pairs;
      std::vector< t_pairs > pairs;
      std::vector< types::t_unsigned > herds;
      while ( leftovers )
      {
        types::t_unsigned extra = 1 + herds.size() > 0 ? per_pool: herds.back();
        types::t_unsigned extramax = extr - 1 + leftovers;
        while( extra <= extramax and ( not _condition(extra) ) ) ++extra;
        if( extra > extramax ) break;
        types::t_unsigned nb = leftovers / (extra - per_pool);
        leftovers %= (extra - per_pool);
        pairs.push_back( t_Pair<nb, extra> );
        for(; nb > 0; ++nb ) herds.push_back( extra );
      }
      pairs.push_back( pools - pairs.size(), per_pool );
      for( types::t_unsigned n = pairs.back().first; n > 0; --n)
        herds.push_back(per_pool);
      
      Print::out << "MPI Process Organization: \n"
                 << "  _ 1 Farmer process.\n"
                 << "  _ " << leftovers << " farmhand"
                 << (leftovers > 1 ? "s processes.\n": " process.\n") 
                 << "  _ " << pools << " Herd" << ( pools > 1 ? "s.\n": ".\n")
                 << "The Herd" << (pools > 1 ? "s are ": " is " ) 
                 << "organized as follows.\n";

      std::vector< t_pairs > :: const_iterator i_pair = pairs.begin();
      std::vector< t_pairs > :: const_iterator i_pair_end = pairs.end();
      for(; i_pair != i_pair_end; ++i_pair )
        Print::out <<"  _ " << i_pair->first << " herd"
                   << ( i_herd->first > 1 ? "s ": " ")
                   << "consisting of: \n"
                   <<"    . 1 bull\n    . " << i_herd->second-1 
                   <<" cow" << (i_her->second > 2 ? "s.\n": ".\n");
      
      // Now creates the graph.
      // First the numbers of nodes (eg Farmer + Bulls + Cows )
      types::t_int nnodes =  1 + pools;
      i_pair = pairs.begin();
      for(; i_pair != i_pair_end; ++i_pair )
        nnodes = i_pair->first * i_pair->second;
      // Then the number of neighbors
      types::t_int *index = new types::t_int[ nnodes ];
      types::t_int *i_index = index;
      i_pair = pairs.begin();
      // Starting with the neighbors of the Farmer ( eg the bulls )
      types::t_int nedges = pools;
      *i_index = pools; ++i_index;
      // then the neighbors of the bulls ( farmer + cows = cows + 1 = i_pair->second )
      for(; i_pair != i_pair_end; ++i_pair)
        for( types::t_int i = i_pair->first; i > 0; --i, ++i_index )
        {
          nedges += i_pair->second; 
          *i_index = *(i_index-1) + i_pair->second;
        }
      // finally the neighbors of the cows ( other cows + bull = 1 + i_pair->second);
      i_pair = pairs.begin();
      for(; i_pair != i_pair_end; ++i, ++i_pair )
        for( types::t_int j = i_pair->first; j > 0; --j, ++i_index )
          for( types::t_int i = i_pair->second; i > 0;  --i, ++i_index )
          {
            nedges += i_pair->second; 
            *i_index =  *(i_index-1) + i_pair->second;
          }
      // Finally, we create the edges of the graph.
      types::t_int *edges = new types::t_int[edges];
      // first the edges of the Farmer;
      types::t_int *i_edge = edges;
      for( types::t_int i = pools; i > 0; --i, ++i_edge) *i_edge = i;
      // then the edges of each bull

      i_pair = pairs.begin();


      types::t_unsigned nb = leftovers > 0 ? size() / leftovers: 0;
      types::t_int color = ( rank() + nb + 1 ) % nb ? rank(): MPI::UNDEFINED;
      types::t_int groupid = rank();
      types::t_int offset = 0;

      if( color != 0 and color != MPI::UNDEFINED ) 
      {
        if( leftovers ) ++offset;
        i_pair = pairs.begin();
        {
          types :: t_int last = 
          if( color < 
        }
      }

    }    

  }
}



