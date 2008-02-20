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
          t_CommBase::wait_bulls();
        
          _offspring.resize(target);   // you might have generated a few more
        }
        template<class T_GATRAITS>
        inline void Farmer<T_GATRAITS> :: onWait( types::t_int _bull )
        {
          types::t_int buff;
          t_Individual indiv;
          t_CommBase::receive_individual( _bull, indiv );
          offspring->push_back( indiv );
          t_CommBase::send_command( _bull, offspring->size() >= target ?
                                             t_CommBase::t_Commands::DONE: 
                                             t_CommBase::t_Commands::GO )
          if( offspring->size() < target ) t_CommBase::activate(_bull);
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
          t_CommBase :: command( t_CommBase::t_CowCommands::DONE );
      
          _offspring.clear();
        }
      } // namespace Breeder
    } // namespace Graph

  } // namespace mpi
} // namespace GA
