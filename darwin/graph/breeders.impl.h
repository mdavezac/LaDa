//
//  Version: $Id$
//

namespace GA
{
  namespace mpi 
  {

    namespace Graph
    {
      namespace Breeders
      {

        template<class T_GATRAITS>
        inline void Farmer<T_GATRAITS> :: operator()(const t_Population& _parents,
                                                     t_Population& _offspring)
        {
          target = (*t_Base::howMany)( (types::t_unsigned) _parents.size());
      
          _offspring.clear();
          offspring = &_offspring;
          t_CommBase::start_all();
        
          while (_offspring.size() < target)
            t_CommBase::test_bulls();
          t_CommBase::wait_bulls();
        
          _offspring.resize(target);   // you might have generated a too many
        }
        template<class T_GATRAITS>
        inline void Farmer<T_GATRAITS> :: onWait( types::t_int _bull )
        {
          types::t_int buff;
          t_Individual indiv;
          t_CommBase :: comm->recv( _bull, ONWAIT_TAG( TAG ), indiv );
          offspring->push_back( indiv );
          t_CommBase::send_command( _bull, offspring->size() >= target ?
                                             t_CommBase::t_Commands::DONE: 
                                             t_CommBase::t_Commands::GO );
          if( offspring->size() < target ) t_CommBase::activate(_bull);
        }
      
        template<class T_GATRAITS>
        inline void Bull<T_GATRAITS> :: operator()(const t_Population& _parents,
                                                   t_Population& _offspring)
        {
          __ASSERT( not t_Base::select, "Selection pointer is not assigned." )
          __ASSERT( not t_Base::age, "Age pointer is not assigned." )
          __ASSERT( not t_Base::op, "Op pointer is not assigned." )
          eoSelectivePopulator<t_Individual> it(_parents, _offspring, *t_Base::select);
          do
          {
            if( not t_Base :: op ) Print :: out << "operator not set" << Print ::endl;
            (*t_Base::op)(it);
            (*it).set_age( (*t_Base::age)() );
            t_CommBase :: request( t_CommBase::t_Requests::WAITING );
            t_CommBase :: comm->send( 0, ONWAIT_TAG( TAG ), *it );
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
