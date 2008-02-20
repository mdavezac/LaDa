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

    namespace Graph
    {
      namespace Breeder
      {

        template<class T_GATRAITS>
        inline void Farmer<T_GATRAITS> :: operator()(const t_Population& _parents,
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
        inline void Farmer<T_GATRAITS> :: onWait( types::t_int _bull )
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
        inline void Farmer<T_GATRAITS> :: onTaboo( types::t_int _bull )
        {
          __ASSERT( taboo, "Taboo pointer has not been set.\n")
          t_Individual indiv;
          t_CommBase::receive_individual( _bull, indiv );
          t_CommBase::send_command( _bull, (*taboo)( indiv ) );
        }
        template<class T_GATRAITS>
        inline void Farmer<T_GATRAITS> :: onObjective( types::t_int _bull )
        {
          __ASSERT( objective, "Objective pointer not set.\n" )
          MPI::DOUBLE buff = 0;
          typename t_Individual :: t_Quantities quantities;
          t_CommBase::receive_quantities( quantities );
          typename t_Individual :: t_Fitness fitness = (*objective)( quantities );
          t_CommBase::send_fitness( fitness );
        }
        template<class T_GATRAITS>
        inline void Farmer<T_GATRAITS> :: onHistory( types::t_int _bull )
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
        inline void Bull<T_GATRAITS> :: operator()(const t_Population& _parents,
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
      } // namespace Breeder
    } // namespace Graph

  } // namespace mpi
} // namespace GA
