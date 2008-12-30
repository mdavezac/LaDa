//
//  Version: $Id$
//

#ifdef _MPI
#include <boost/mpi/collectives.hpp>
#endif

namespace LaDa
{
  namespace GA
  {
    
 
    template< class T_STORE, class T_EVALUATION >
    inline void PrintGA<T_STORE,T_EVALUATION> :: operator()()
    {
      if ( not ( do_print_each_call or store.newresults() ) )
      {
        Print::xmg << Print::Xmg::clearall;
        return;
      }
 
      std::string special = "";
      if ( not store.newresults() ) special = " ? ";
     
      Print::xmg << Print::Xmg::comment << special << "Iteration " << age() 
                 << Print::endl 
                 << Print::Xmg::comment << special << "Evaluation Calls: " 
                 << evaluation.nb_eval << " " << evaluation.nb_grad 
                 << Print::endl;
      
      if( store.newresults() ) store.print_results( age() );
      Print::xmg << Print::flush;
    }
 
    template< class T_STORE, class T_EVALUATION >
    inline void PrintGA<T_STORE,T_EVALUATION> :: lastCall()
    {
      Print::xmg << Print::Xmg::comment << "Last Found Result" << Print::endl;
      store.print_results(age(), true);
      Print::xmg << Print::flush;
    }
 
 
 
 
    template< class T_GATRAITS>
    inline void PrintOffspring<T_GATRAITS> :: operator()( const t_Population &_pop )
    {
      types::t_unsigned ga_age = age();
      if ( ga_age == 0 ) return;
      typename t_Population :: const_iterator i_pop = _pop.begin();
      typename t_Population :: const_iterator i_end = _pop.end();
 
      // The population  of a sorted stat is not an eoPop, but rather a vector of
      // pointers to individuals in an eoPop. Hence the double deferencing
      Print::out << "New Individuals: \n";
      for( ; i_pop != i_end; ++i_pop )
        if ( ga_age == (*i_pop)->get_age() )
          Print::out << *(*i_pop)  << "  Fitness: " 
                     << Print::fixed << Print::setw(12) << Print::setprecision(5) << "  "
                     <<  (*i_pop)->fitness() << "\n";
      Print::out << "\n";
    }
 
 
 
    template< class T_GATRAITS>
    inline void PrintPop<T_GATRAITS> :: operator()( const t_Population &_pop )
    {
      typename t_Population :: const_iterator i_pop = _pop.begin();
      typename t_Population :: const_iterator i_end = _pop.end();
 
      // The population  of a sorted stat is not an eoPop, but rather a vector of
      // pointers to individuals in an eoPop. Hence the double deferencing
      Print::out << "Current Population: \n";
      for( ; i_pop != i_end; ++i_pop )
        Print::out << *(*i_pop)  << "  Fitness: " 
                   << Print::fixed << Print::setw(12) << Print::setprecision(5) << "  "
                   <<  (*i_pop)->fitness() << "\n";
      Print::out << "\n";
    }
 
 
    template< class T_VALUE, class T_BINOP, class T_GATRAITS>
    inline bool Terminator<T_VALUE, T_BINOP, T_GATRAITS> :: operator() (const t_Population &_pop )
    {
      if ( (not term ) or ( term and binop( ref, term ) ) )
       return true;
      lastCall();
      return false;
    }
 
    template< class T_VALUE, class T_BINOP, class T_GATRAITS>
    inline void Terminator<T_VALUE, T_BINOP, T_GATRAITS> :: lastCall() const
    {
      Print::xmg << Print::Xmg::comment << "Terminator, type: " << type
                 << ", ref= " << ref 
                 << ", term=" << term << Print::endl;
      Print::out << "\n Stopping run on:\n    Terminator, type: " << type
                 << ", ref= " << ref 
                 << ", term=" << term << Print::endl;
    }
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
    template< class T_GATRAITS>
    inline void IslandsContinuator<T_GATRAITS> :: apply_stats(const t_Population& _pop) const
    {
      // first goes through sorted stats
      if ( not sorted.empty() )
      { 
        typename std::list < eoSortedStatBase<t_Individual>* > :: const_iterator
             i_sorted = sorted.begin();
        typename std::list < eoSortedStatBase<t_Individual>* > :: const_iterator
             i_end = sorted.end();
        std::vector<const t_Individual*> sorted_pop;
        _pop.sort(sorted_pop);
        for ( ; i_sorted != i_end; ++i_sorted )
          (*i_sorted)->operator()(sorted_pop);
      }
      // then through unsorted stats
      if ( not stats.empty() )
      { 
        typename std::list < eoStatBase<t_Individual>* > :: const_iterator
               i_stat = stats.begin();
        typename std::list < eoStatBase<t_Individual>* > :: const_iterator
               i_end = stats.end();
        for ( ; i_stat != i_end; ++i_stat )
          (*i_stat)->operator()(_pop);
      }
    } 
 
    template< class T_GATRAITS>
    inline void IslandsContinuator<T_GATRAITS> :: apply_monitors_updaters()
    {
      // applies monitors
      if ( not monitors.empty() )
      { 
        typename std::list < eoMonitor* > :: const_iterator i_monitor = monitors.begin();
        typename std::list < eoMonitor* > :: const_iterator i_end = monitors.end();
        for ( ; i_monitor != i_end; ++i_monitor )
          (*i_monitor)->operator()();
      }
      // then through updaters
      if ( not updaters.empty() )
      { 
        typename std::list < eoUpdater* > :: const_iterator i_updater = updaters.begin();
        typename std::list < eoUpdater* > :: const_iterator i_end = updaters.end();
        for ( ; i_updater != i_end; ++i_updater )
          (*i_updater)->operator()();
      }
    }
 
    template< class T_GATRAITS>
    inline bool IslandsContinuator<T_GATRAITS>
      :: apply_continuators( const t_Population & _pop ) const
    {
      typename std::list < eoContinue<t_Individual>* > :: const_iterator
           i_continue = continuators.begin();
      typename std::list < eoContinue<t_Individual>* > :: const_iterator
           i_end = continuators.end();
      bool result = true;
      for (; i_continue != i_end; ++i_continue )
        if ( not (*i_continue)->operator()( _pop) )
          result = false;
      return result;
    }
 
    template< class T_GATRAITS>
    inline void IslandsContinuator<T_GATRAITS> :: lastCall( const t_Population& _pop ) 
    {
      // first goes through sorted stats
      if ( not sorted.empty() )
      { 
        typename std::list < eoSortedStatBase<t_Individual>* > :: const_iterator 
            i_sorted = sorted.begin();
        typename std::list < eoSortedStatBase<t_Individual>* > :: const_iterator
            i_end = sorted.end();
        std::vector<const t_Individual*> sorted_pop;
        _pop.sort(sorted_pop);
        for ( ; i_sorted != i_end; ++i_sorted )
          (*i_sorted)->lastCall(sorted_pop);
      }
      // then through unsorted stats
      if ( not stats.empty() )
      { 
        typename std::list < eoStatBase<t_Individual>* > :: const_iterator i_stat = stats.begin();
        typename std::list < eoStatBase<t_Individual>* > :: const_iterator i_end = stats.end();
        for ( ; i_stat != i_end; ++i_stat )
          (*i_stat)->lastCall(_pop);
      }
      // applies monitors
      if ( not monitors.empty() )
      { 
        typename std::list < eoMonitor* > :: const_iterator i_monitor = monitors.begin();
        typename std::list < eoMonitor* > :: const_iterator i_end = monitors.end();
        for ( ; i_monitor != i_end; ++i_monitor )
          (*i_monitor)->lastCall();
      }
      // then through updaters
      if ( not updaters.empty() )
      { 
        typename std::list < eoUpdater* > :: const_iterator i_updater = updaters.begin();
        typename std::list < eoUpdater* > :: const_iterator i_end = updaters.end();
        for ( ; i_updater != i_end; ++i_updater )
          (*i_updater)->lastCall();
      }
    }
 
    template< class T_GATRAITS>
    inline bool IslandsContinuator<T_GATRAITS> :: apply ( iterator &_i_begin, iterator &_i_end )
    {
      iterator i_pop = _i_begin;
 
      for( ; i_pop != _i_end; ++i_pop )
        apply_stats( *i_pop );
 
      apply_monitors_updaters();
 
      Print::xmg << Print::flush; 
 
      bool result =    ( max_generations and generation_counter() < max_generations ) 
                    or ( not max_generations );
      for( i_pop = _i_begin; i_pop != _i_end; ++i_pop )
        if ( not apply_continuators( *i_pop ) )
          result = false;
 
          // checks if stop file exists
      __MPICODE
      (
        boost::mpi::communicator world;
        if ( world.rank() == 0 ) 
      )
        {
          std::ifstream file( stop_filename.c_str(), std::ios_base::in );
          if ( file.is_open() )
          {
            Print::xmg << Print::Xmg::comment << "Stopping on finding file "
                       << stop_filename << Print::endl;
            Print::out << "\n\nStopping on finding file "
                       << stop_filename << Print::endl;
            result = false; 
          }
        }
      __MPICODE( result = boost::mpi::all_reduce( world, result, std::logical_and<bool>() ); )
 
      // last call
      if ( not result )
      {
        for( i_pop = _i_begin; i_pop != _i_end; ++i_pop )
          lastCall( *i_pop );
      }
 
      // increments generation
      ++generation_counter;
 
      return result;
    }
 
 
 
 
 
 
    template < class T_CLASS >
    inline void SaveEvery<T_CLASS> :: operator()()
    {
      ++age;
 
      if ( age % every )
        return;
 
      if ( (object.*func)() )
        return;
 
      std::cerr << "Could not perform save " << std::endl;
    }
 
 
 
 
 
 
 
#  ifdef _MPI
    template< class T_TYPE >
    inline void Synchronize<T_TYPE> :: operator()()
    {
      t_Type diff = object - current_value;
      boost::mpi::communicator world;
      diff = boost::mpi::all_reduce( world, diff, std::plus<t_Type>() );
      current_value += diff;
      object = current_value;
    }
#  endif

} // namespace LaDa

}
