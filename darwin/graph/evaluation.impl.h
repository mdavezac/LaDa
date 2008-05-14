//
//  Version: $Id$
//



namespace GA
{
  namespace mpi 
  {
    namespace Graph
    {
      namespace Evaluation
      {
        template<class T_GATRAITS, template< class > class T_BASE >
          void Farmer<T_GATRAITS, T_BASE> :: operator()( t_Population &_pop,
                                                         t_Population &_offspring ) 
            { 
              __ASSERT(     t_Base :: objective
                        and t_CommBase :: objective, 
                        "objective pointer has not been set.\n" )
              // eventually, this will call Evaluator::Farmer::evaluate(), 
              // and put unknown individuals into a list.
                Print::out << "evaluation: 0" << Print::endl;
              t_Base::evaluate(_offspring); 

                Print::out << "evaluation: 1" << Print::endl;
                cache_eval.print();

              // At this point,  unknowns are evaluated. 
              t_CommBase :: start_all();
              while( cache_eval.notdone() ) t_CommBase::test_bulls();
              t_CommBase :: wait_bulls();

              // Continue with real evaluations.
              if ( t_Base::objective->is_valid() ) return;
              
              // if invalid, recomputes whole population
              t_Base::evaluate( _pop );
              t_Base::evaluate( _offspring );
            }

        template<class T_GATRAITS, template< class > class T_BASE >
          void Farmer<T_GATRAITS, T_BASE> :: onWait( types::t_unsigned _bull )
          {
            Print::out << "onWait" << Print::endl;
            t_Individual *next = cache_eval.onWait( _bull );
            if( not next )
            {
              t_CommBase::send_command( _bull, t_CommBase::t_Commands::DONE );
              return;
            }
            Print::out << "onWait 1 " << Print::endl;
            t_CommBase::send_command( _bull, t_CommBase::t_Commands::GO );
            Print::out << "onWait 2 " << Print::endl;
            t_CommBase::send_individual( _bull, *next );
            Print::out << "onWait 3 " << Print::endl;
            t_CommBase::activate( _bull );
            Print::out << "onWait 4 " << Print::endl;
          }

        template<class T_GATRAITS, template< class > class T_BASE >
          void Bull<T_GATRAITS, T_BASE> :: operator()(t_Population& _parents,
                                                      t_Population& _offspring)
          {
            Print::out << "evaluate( t_pop, t_pop)" << Print::endl;
              t_CommBase :: request( t_CommBase::t_Requests::WAITING );
            Print::out << "evaluate( t_pop, t_pop) -1" << Print::endl;
            while ( t_CommBase::obey() != t_CommBase :: t_Commands :: DONE )
            {
            Print::out << "evaluate( t_pop, t_pop) obey 0 " << Print::endl;
              t_Individual individual;
              t_CommBase :: receive_individual( individual );
            Print::out << "evaluate( t_pop, t_pop) obey 1 " << Print::endl;
              metaeval.init( individual );
            Print::out << "evaluate( t_pop, t_pop) obey 2 " << Print::endl;
              metaeval.evaluate();
            Print::out << "evaluate( t_pop, t_pop) obey 3 " << Print::endl;
              t_CommBase :: request( t_CommBase::t_Requests::WAITING );
            Print::out << "evaluate( t_pop, t_pop) obey 4 " << Print::endl;
            }
            Print::out << "evaluate( t_pop, t_pop) end " << Print::endl;
          }

      } // namespace Evaluation
    } // namespace Graph
  } // namespace mpi
} // namespace GA
