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
          void Farmer<T_GATRAITS, T_BASE>
            :: operator()( t_Population &_pop, t_Population &_offspring ) 
            { 
              __ASSERT( objective, "objective pointer has not been set.\n" )
              // eventually, this will call Evaluator::Farmer::evaluate(), 
              // and put unknown individuals into a list.
              t_Base::evaluate(_offspring); 
            
              // At this point,  unknowns are evaluated. 
              t_CommBase :: start_all();
              while( unknowns.size() ) t_CommBase::test_bulls();
              t_CommBase :: wait_bulls();

              // Continue with real evaluations.
              if ( objective->is_valid() ) return;
              
              // if invalid, recomputes whole population
              t_Base::evaluate( _pop );
              t_Base::evaluate( _offspring );
            }

        template<class T_BASE>
        void Farmer<T_BASE> :: onWait( types::t_unsigned _bull )
        {
          t_Individual *next = cache_eval.onWait( _bull );
          if( not next )
          {
            t_CommBase::send_command( _bull, t_CommBase::t_Commands::DONE );
            return;
          }
          t_CommBase::send_command( _bull, t_CommBase::t_Commands::GO );
          t_CommBase::send_individual( _bull, *next );
          i_unknwon->second = _bull;
          t_CommBase::activate( _bull );
        }

        template<class T_BASE>
        inline void Bull<T_BASE> :: operator()(const t_Population& _parents,
                                               t_Population& _offspring)
        {
          while ( t_CommBase::obey() != t_CommBase :: t_Commands :: DONE )
          {
            t_Individual individual;
            t_CommBase :: receive_individual( individual );
            evaluator.init( individual );
            evaluator.evaluate( individual );
            t_CommBase :: request( t_CommBase::t_Requests::WAITING );
          }
        }

      } // namespace Evaluation
    } // namespace Graph
  } // namespace mpi
} // namespace GA
