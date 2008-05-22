//
//  Version: $Id$
//

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

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
      template< class _ITERATOR >
        void print_edges( _ITERATOR _first, _ITERATOR _last )
        {
          for(; _first != _last; ++_first )
            Print::out << _first->first << " is connected to " 
                       << _first->second << "\n"; 
        }

      template< class T_CONDITION >
      bool Topology :: init( T_CONDITION _condition )
      {
        try
        {
          if ( pools < 2 ) return false;
          types::t_int per_pool = ( comm.size() - 1) / pools;
          if ( per_pool <= 0 ) return false;
          
          typedef std::pair<types::t_int, types::t_int > t_Edge;
          typedef std::vector< t_Edge > t_Edges;
          t_Edges edges;
          // first add farmer-bull edges
          for( types::t_int i = pools; i > 0; --i )
            edges.push_back( t_Edge( 0, i ) );
          // Then add for each bull, add the cow edges.
          for( types::t_int j= comm.size(); j > pools; )
            for( types::t_int i = pools; i > 0 and j > pools; --i, --j )
              edges.push_back( t_Edge(i, j) );
          
          // Now creates graph.
          typedef boost::adjacency_list< boost::vecS, 
                                         boost::vecS, 
                                         boost::bidirectionalS > t_Graph;
          t_Graph g( edges.begin(), edges.end(), comm.size() );
          
          // Creates a graph_communicator
          if( graph_comm ) delete graph_comm;
          __TRYCODE( 
            graph_comm = new boost::mpi::graph_communicator( comm, g, true );,
            "Could not create graph communicator.\n" 
          )
          __DOASSERT( not graph_comm, "Could not create graph communicator.\n")
          
          // Finally, prints graph
          Print::out << "Graph Connection:\n";
          print_edges( boost::mpi::edges(*graph_comm).first, 
                       boost::mpi::edges(*graph_comm).second );
          Print::out << Print::endl;
          
          
          // Now creates groups
          types::t_int color = MPI::UNDEFINED;
          if( ((::MPI::Graphcomm) *graph_comm) == MPI::COMM_NULL )
            { type = t_Type::t_Type::FARMHAND; goto exit; }
          
          // communicator between Farmer and Bulls.
          if( comm.rank() == 0 ) type = t_Type::FARMER;
          if( comm.rank() <= pools ) color = 0;
          head_comm = comm.split( color );
          // Now we create the communicator for the herd.
          color = MPI::UNDEFINED;
          if( comm.rank() <= pools and comm.rank() > 0 )
          {
            color = comm.rank();
            type = t_Type::BULL;
          }
          for( types::t_int j= comm.size(); j > pools; )
            for( types::t_int i = pools; i > 0 and j > pools; --i, --j )
              if( comm.rank() == j )
              {
                color = i;
                type = t_Type::COW;
              }
          pool_comm = comm.split( color );
        }
        catch( std::exception &_e )
        {
          std::ostringstream sstr;
          sstr << __SPOT_ERROR << "Will not create graph topology:\n"
               << _e.what() << "\n";
          std::cerr << sstr.str();
          Print :: out << sstr.str();
          return false;
        }
 
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
          sstr << "." << pool_comm.rank();
          _eval.set_mpi( &pool_comm, sstr.str() );
        }

    } // namspace Graph
      
  } // namspace mpi

} // namespace GA
