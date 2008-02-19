//
//  Version: $Id$
//
#ifndef  _DARWIN_COMMUNICATORS_H_
#define  _DARWIN_COMMUNICATORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace GA
{
  namespace mpi 
  {

    //! \brief Creates a GA mpi topology consisting of a single unit.
    template< class T_CONDITION >
    void Graph<T_CONDITION>::init_collectivists( T_CONDITION &_condition )
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

    //! \brief Creates the GA mpi topology and groups.
    //! \details The topology and subsequent groups can be explained either
    //           poetically or mathematically. One could imagine a farm headed
    //           by a lonely single man. He believes himself intrusted by God
    //           with the propagation of the \e Bos \e Taurus. He was granted a
    //           number of bulls. And each bull a number of cows with which to
    //           do the will of the Almighty. In his wisdom, the Master
    //           Universe may have afflicted good-for-nothing farmhands to the
    //           good man. More specifically, if there are \f$ n = p*h + 1\f$
    //           processors. Processor +1 is the head boss directing \e p pools
    //           of \e h processor. The processors of a pool work together to
    //           evaluate a configuration. In any case, the topology consists
    //           of an origin (the farmer) and \p rings. The origin is
    //           connected to a single member of each ring. That's the bull.
    //           The cows make up the other members of the ring. Rings
    //           consisting of asingle bull are allowed (ring becomes a point),
    //           as well as rings consisting of a cow and a bull ( ring becomes
    //           a segment). 
    //
    //           There are three "inputs" to the graph: \e n + 1 the number of
    //           processes, \p the requested number of pools, \a _condition a
    //           condition which \e h must fulfill (for instance, \t escan
    //           requires that the number of waves is divisible by the number of
    //           procs in a pool). Hence the farmhands. The algorithm strives
    //           to find the best combination such that \b at \b most \e p pools
    //           are created, and that these pools have a balanced number of
    //           processor. It applies the condition (\a _condition) \e sine \e
    //           qua \e non.
    template< class T_CONDITION >
    void Graph<T_CONDITION>::init( T_CONDITION &_condition )
    {
      //! One single pool
      if( pools <= 1 ) init_collectivists( _condition );
      types::t_int per_pool = ( size() - 1) / pools;
      types::t_int leftovers = ( size() - 1) % pools;

      while ( per_pool and ( not _condition( per_pool) ) ) --per_pool;
      leftovers += size() - per_pool * pools;

      typedef std::pair<types::t_int, types::t_int > t_pairs;
      std::vector< types::t_pairs > pairs(1, t_pairs(pools, per_pool) );
      std::vector< types::t_int > herds( pools, per_pool);
      while ( leftovers )
      {
        types::t_unsigned extra = 1 + herds.front();
        types::t_unsigned extramax = herds.front() + leftovers;
        while( extra < extramax and ( not _condition(extra) ) ) ++extra;
        types::t_unsigned nb = leftovers / (extra - per_pool);
        if( extra >= 2 * herds.front() and 2 * nb < pools )
        { 
          std::cerr << "Herds cannot be balanced.\n"
                    << "reducing the number of herds.\n";
          --pool;
          init( _condition );
          return;
        }
        if( extra == extramax ) break;
        leftovers %= (extra - per_pool);
        std::fill(herds.begin(), herds.begin() + nb, extra);
        pairs.front().first -= nb;
        pairs.push_back( t_pairs( nb, extra ) );
      }
      
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
                   << ( i_pair->first > 1 ? "s ": " ")
                   << "consisting of: \n"
                   <<"    . 1 bull\n    . " << i_pair->second-1 
                   <<" cow" << (i_her->second > 2 ? "s.\n": ".\n");
      
      // Now creates the graph.
      // First the numbers of nodes (eg Farmer + Bulls + Cows )
      types::t_int nnodes =  1;
      i_pair = pairs.begin();
      for(; i_pair != i_pair_end; ++i_pair )
        nnodes = i_pair->first * i_pair->second;
      // Then the number of neighbors.
      types::t_int *indices = new types::t_int[ nnodes ];
      types::t_int *i_index = indices;
      types::t_int *i_index_end = indices + nnodes;
      // Starting with the neighbors of the Farmer, connected to each bull
      *i_index = pools; ++i_index;
      // then the neighbors of the bulls, connected to the farmer, and to two
      // cows. The complete herd forms a ring
      // EXCEPTION CASES: 1. There are no cows in the herd.
      //                  2. There is only one cow in the herd.
      std::vector<types::t_unsigned> :: const_iterator i_herd = herd.begin();
      std::vector<types::t_unsigned> :: const_iterator i_herd_end = herd.end();
      for(; i_herd != i_herd_end; ++i_herd, ++i_index )
        switch( *i_herd )
        {
          case 1:  *i_index = 1 + *(i_index - 1); break;
          case 2:  *i_index = 2 + *(i_index - 1); break;
          default: *i_index = 3 + *(i_index - 1); break;
        }
      // Each cow is connected to two other members of its herd.
      // EXCEPTION CASES: 1. There are no cows in the herd.
      //                  2. There is only one cow in the herd.
      // The first exception case does not enter the second loop.
      i_herd = herd.begin();
      for(; i_herd != i_herd_end; ++i_herd )
        for( types::t_int i = *i_herd - 1; i > 0; --i, ++i_index )
          *i_index = *(i_index-1) + ( *i_herd == 2 ? 1: 2); 
      // Finally, we create the edges of the graph.
      types::t_int *edges = new types::t_int[ *(i_index - 1) ];
      // first the edges of the Farmer.
      types::t_int *i_edge = edges;
      for( types::t_int i = pools; i > 0; --i, ++i_edge) *i_edge = i;
      // then the edges of each bull.
      i_herd = herds.begin();
      types::t_int first_cow_in_herd = pools + 1;
      for(; i_herd != i_herd_end; ++i_herd )
      {
        *i_edge = 0; ++i_edge;
        if( *i_herd == 1 ) continue;
        *i_edge = first_cow_in_herd; ++i_edge;
        if( *i_herd == 2 ) continue;
        *i_edge = first_cow_in_herd + (*i_herd) - 2; ++i_edge;
      }
      // now for the cows
      // First cow is linked to the bull and the second cow.
      // Last cow is linkerd to the previous cow and the bull.
      // All other cows are linked to the next and previous cow.
      i_herd = herd.begin();
      types::t_int current_cow = pools + 1;
      for(types::t_int i = 1; i_herd != i_herd_end; ++i_herd, ++i )
      {
        if( *i_herd == 1 ) continue;
        *i_edge = i; ++i_edge; 
        if( *i_herd == 2 ) continue;
        *i_edge = current_cow + 1; ++i_edge; ++current_cow;
        for( types::t_int i = *i_herd - 3; i > 0; --i, ++next_cow)
        {
          *i_edge = current_cow - 1; ++i_edge;
          *i_edge = current_cow + 1; ++i_edge;
        }
        *i_edge = current_cow - 1; ++i_edge; ++current_cow;
        *i_edge = i; ++i_edge; 
      }


      // At this point, we can construct the graph
      graph_comm =  comm->Create_graph( nnodes, indices, edges, true);
      if( graph_comm = MPI::COMM_NULL ) type = FARMHAND;

      // Now we can create the herds.
      // First we create the communicator between Farmer and Bulls.
      types::t_int color = MPI::UNDEFINED;
      if( rank() == 0 ) type = FARMER;
      if( rank() < pools ) color = 0;
      farmer_comm = &comm->Split( color, rank() );
      // Now we create the communicator for the herd.
      color = MPI::UNDEFINED;
      if( rank() < pools and rank() > 0 ) { color = rank(); type = BULLS; }
      first_cow_in_herd = pools + 1;
      for(types::t_int i = 1; i_herd != i_herd_end; ++i_herd, ++i )
        if(     rank() >= first_cow_in_herd
            and rank() < first_cow_in_herd - (*i_herd) - 1 )
        {
          color = i;
          type = COW;
        }
      herd_comm = &comm->Split( color, rank() );

    }    
  }


  template<class T_GATRAITS>
  inline Graph::Breeder<T_GATRAITS>* Graph :: create_breeder( eoState &_state )
  {
    Breeder * result = NULL;

    try
    {
      switch( topo->type )
      {
        case mpi::Graph::COLLECTIVISTS: 
           result = new Breeder::Collectivist( this ); break;
        case mpi::Graph::COW: 
           result = new Breeder::Cow( this ); break;
        case mpi::Graph::BULL: 
           result = new Breeder::Bull( this ); break;
        case mpi::Graph::FARMHANDS: return; break;
           result = new Breeder::FarmhandGraphBreeder( this ); break;
      }

      _state.saveFunctor( result );
      return result;
    }
    catch ( std::exception &_e ) 
    {
      if( result ) delete result;
      __THROW_ERROR( "Error encountered while creating Breeder.\n" << _e.what() )
    }
  }
  

  template<class T_GATRAITS>
  inline void FarmerGraphBreeder<T_GATRAITS> 
    :: operator()(const t_Population& _parents,
                  t_Population& _offspring)
  {
    target = (*howMany)( (types::t_unsigned) _parents.size());

    _offspring.clear();
    offspring = &_offspring;
    t_CommBase::startall();
  
    while (_offspring.size() < target)
      t_CommBase::test_bulls();
  
    _offspring.resize(target);   // you might have generated a few more
  }
  template<class T_GATRAITS>
  inline void FarmerGraphBreeder<T_GATRAITS> :: onWait( types::t_int _bull )
  {
    types::t_int buff;
    t_Individual indiv;
    t_CommBase::receive_individual( _bull, indiv );
    offspring->push_back( indiv );
    t_CommBase::send_command( offspring->size() >= target ? t_CommBase::DONE: 
                                                            t_CommBase::Done, 
                               _bull );
  }
  template<class T_GATRAITS>
  inline void FarmerGraphBreeder<T_GATRAITS> :: onTaboo( types::t_int _bull )
  {
    __ASSERT( taboo, "Taboo pointer has not been set.\n")
    t_Individual indiv;
    t_CommBase::receive_individual( _bull, indiv );
    t_CommBase::send_command( _bull, (*taboo)( indiv ) );
  }
  template<class T_GATRAITS>
  inline void FarmerGraphBreeder<T_GATRAITS> :: onObjective( types::t_int _bull )
  {
    __ASSERT( objective, "Objective pointer not set.\n" )
    MPI::DOUBLE buff = 0;
    typename t_Individual :: t_Quantities quantities;
    t_CommBase::receive_quantities( quantities );
    typename t_Individual :: t_Fitness fitness = (*objective)( quantities );
    t_CommBase::send_fitness( fitness );
  }
  template<class T_GATRAITS>
  inline void FarmerGraphBreeder<T_GATRAITS> :: onHistory( types::t_int _bull )
  {
    // Don't expect this message if history is not set
    __ASSERT( history, "History Pointer not set.\n")
    t_Individual indiv;
    t_CommBase::receive_individual( _bull, indiv );
    MPI::BOOL buff = history->clone( indiv );
    t_CommBase::send_command( buff, _bull );
    if( not buff ) return;
    t_CommBase::send_individual( _bull, indiv );
  }



  template<class T_GATRAITS>
  inline void BullGraphBreeder<T_GATRAITS> 
    :: operator()(const t_Population& _parents,
                  t_Population& _offspring)
  {
    eoSelectivePopulator<t_Individual> it(_parents, _offspring, select);
    do
    {
      (*op)(it);
      (*it).set_age( age() );
      t_CommBase :: request( t_CommBase::t_Requests::WAITING );
      t_CommBase :: send_individual();
      ++it;
    }
    while ( t_CommBase::obey() != t_CommBase :: t_Commands :: DONE );
    t_CommBase :: bcast( t_CommBase::t_CowCommands::DONE );

    _offspring.clear();
  }

}
