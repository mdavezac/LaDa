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
        template<class T_BASE>
        void Farmer<T_BASE> :: operator()(const t_Population& _parents,
                                          t_Population& _offspring)
        {
          __ASSERT( evaluator, "Pointer to Evaluator is not set.\n")
          __ASSERT( objective, "Pointer to Evaluator is not set.\n")
          __ASSERT( store, "Pointer to Evaluator is not set.\n")

          unknowns.clear();
          typename t_Population :: iterator i_indiv = _offspring.begin();
          typename t_Population :: iterator i_indiv_end = _offspring.end();
          for(; i_indiv != i_indiv_end; ++i_indiv )
            if( i_indiv->invalid() ) unknowns.push_back( t_Unknown( i_indiv, -1 ) );

          t_CommBase :: start_all();
          while( unknowns.size() ) t_CommBase::test_bulls();
          t_CommBase :: wait_bulls();

          // Sometimes, evaluation is a bit more complex.
          // Say for instance in case of moving target.
          T_BASE::operator()( _parents, _offpsring );
        }

        template<class T_BASE>
        void Farmer<T_BASE> :: onWait( types::t_unsigned _bull )
        {
          typename t_Unknowns :: iterator i_unknown = unknowns.begin();
          typename t_Unknowns :: iterator i_unknown_end = unknowns.end();
          for(; i_unknown != i_unknown_end; ++i_unknown )
            if( i_unknown->second == _bull - 1 )
            {
              unkowns.erase( i_unknown );
              break;
            }
          i_unknown = unknowns.begin();
          for(; i_unknown != i_unknown_end; ++i_unknown )
            if( i_unknown->second == -1 ) break;

          if( i_unknown == i_unknown_end )
          {
            t_CommBase::send_command( _bull, t_CommBase::t_Commands::DONE );
            return;
          }
          t_CommBase::send_command( _bull, t_CommBase::t_Commands::GO );
          t_CommBase::send_individual( _bull, *i_unknwon->first );
          i_unknwon->second = _bull;
          t_CommBase::activate( _bull );
        }

        //! Sets history for t_Base == Evaluation::WithHistory.
        template<class T_BASE> template <>
        inline void Farmer<T_BASE> 
          :: set_base_history<Evaluation::WithHistory>( t_History *_history )
            { t_Base :: history = _history; }
          

        template<class T_BASE>
        inline void Bull<T_BASE> :: operator()(const t_Population& _parents,
                                               t_Population& _offspring)
        {
          while ( t_CommBase::obey() != t_CommBase :: t_Commands :: DONE )
          {
            t_Individual individual;
            t_CommBase :: receive_individual( individual );
            t_CommBase :: bcast( t_CommBase::t_CowCommands::GO );
            t_CommBase :: bcast( individual );
            evaluate( individual );
            t_CommBase :: request( t_CommBase::t_Requests::WAITING );
          }
        }

        template<class T_BASE>
        inline void Bull<T_BASE> :: evaluate( t_Individual &_indiv)
        {
          while ( t_CommBase::obey() != t_CommBase :: t_Commands :: DONE )
          {
            t_Individual individual;
            t_CommBase :: receive_individual( individual );
            t_CommBase :: bcast( t_CommBase::t_CowCommands::GO );
            t_CommBase :: bcast( individual );
            evaluate( individual );
            t_CommBase :: request( t_CommBase::t_Requests::WAITING );
          }
        }

      } // namespace Evaluation
    } // namespace Graph
  } // namespace mpi
} // namespace GA
