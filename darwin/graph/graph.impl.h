//
//  Version: $Id$
//

#ifdef _DEBUG
#include <boost/lambda/lambda.hpp>
#include <algorithm>
#endif

namespace GA
{
  namespace mpi 
  {

    namespace Graph
    {

      template< class T_CONDITION >
      bool Topology :: init( T_CONDITION _condition )
      {
        if( pools < 2 ) return false;
        types::t_int per_pool = ( size() - 1) / pools;
 
        while ( per_pool and ( not _condition( per_pool) ) ) --per_pool;
        types::t_int leftovers = size() - per_pool * pools - 1;
 
        typedef std::pair<types::t_int, types::t_int > t_pairs;
        std::vector< t_pairs > pairs(1, t_pairs(pools, per_pool) );
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
            --pools;
            return init( _condition );
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
                     <<" cow" << (i_pair->second > 2 ? "s.": ".")
                     << Print::endl;
        
        // Now creates the graph.
        // First the numbers of nodes (eg Farmer + Bulls + Cows )
        types::t_int nnodes =  1;
        i_pair = pairs.begin();
        for(; i_pair != i_pair_end; ++i_pair )
          nnodes += i_pair->first * i_pair->second;
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
        std::vector<types::t_int> :: const_iterator i_herd = herds.begin();
        std::vector<types::t_int> :: const_iterator i_herd_end = herds.end();
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
        i_herd = herds.begin();
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
          first_cow_in_herd += *i_herd - 1;
          if( *i_herd == 2 ) continue;
          *i_edge = first_cow_in_herd - 1; ++i_edge;
        }
        // now for the cows
        // First cow is linked to the bull and the second cow.
        // Last cow is linkerd to the previous cow and the bull.
        // All other cows are linked to the next and previous cow.
        i_herd = herds.begin();
        types::t_int current_cow = pools + 1;
        for(types::t_int i = 1; i_herd != i_herd_end; ++i_herd, ++i )
        {
          if( *i_herd == 1 ) continue;
          *i_edge = i; ++i_edge; 
          if( *i_herd == 2 ) continue;
          *i_edge = current_cow + 1; ++i_edge; ++current_cow;
          for( types::t_int j = *i_herd - 3; j > 0; --j, ++current_cow)
          {
            *i_edge = current_cow - 1; ++i_edge;
            *i_edge = current_cow + 1; ++i_edge;
          }
          *i_edge = current_cow - 1; ++i_edge; ++current_cow;
          *i_edge = i; ++i_edge; 
        }
 
 
        types::t_int color = MPI::UNDEFINED;
        // At this point, we can construct the graph
        graph_comm =  comm->Create_graph( nnodes, indices, edges, true);
        if( graph_comm == MPI::COMM_NULL )
          { type = t_Type::t_Type::FARMHAND; goto exit; }
 
        // Now we can create the herds.
        // First we create the communicator between Farmer and Bulls.
        if( rank() == 0 ) type = t_Type::FARMER;
        if( rank() <= pools ) color = 0;
        head_comm = comm->Split( color, rank() );
        // Now we create the communicator for the herd.
        color = MPI::UNDEFINED;
        if( rank() <= pools and rank() > 0 ) { color = rank(); type = t_Type::BULL; }
        first_cow_in_herd = pools + 1;
        i_herd = herds.begin();
        for(types::t_int i = 1; i_herd != i_herd_end; ++i_herd, ++i )
        {
          if(     rank() >= first_cow_in_herd
              and rank() < first_cow_in_herd + (*i_herd) - 1 )
          {
            color = i;
            type = t_Type::COW;
          }
          first_cow_in_herd += *i_herd - 1;
        }
        pool_comm = comm->Split( color, rank() );
 
 
exit:
#ifdef _DEBUG
         switch( type )
         {
           case t_Type::FARMER:
             Print::out << "This process is the Farmer."; break;
           case t_Type::FARMHAND:
             Print::out << "This process is a Farmhand."; break;
           case t_Type::BULL:
             Print::out << "This process is a Bull."; break;
           case t_Type::COW:
             Print::out << "This process is a Cow."; break;
         }
         Print::out << Print::endl;
#endif
        return true;
      }    

      template <class T_EVALUATOR> 
        void Topology :: set_mpi( T_EVALUATOR &_eval )
        {
          if(    type == t_Type::FARMER 
              or type == t_Type::FARMHAND ) return;
          std::ostringstream sstr;
          ::mpi::Base comm( pool_comm );
          sstr << "." << pool_comm.Get_rank();
          _eval.set_mpi( &comm, sstr.str() );
        }

    } // namspace Graph
      
  } // namspace mpi

} // namespace GA
