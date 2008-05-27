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
              t_Base::evaluate(_offspring); 


              // At this point,  unknowns are evaluated. 
              t_CommBase :: start_all();
              while( notdone() ) t_CommBase::test_bulls();
              Print :: out << "Waiting on all bulls" << Print::endl;
              t_CommBase::comm->barrier();
              Print :: out << "passed barrier" << Print::endl;
//             t_CommBase :: wait_bulls();

              // Ok, now recomputes everything through base class. 
              // if invalid, recomputes whole population
              t_Base::operator()( _pop, _offspring );
            }
        template<class T_GATRAITS, template< class > class T_BASE >
          typename Farmer<T_GATRAITS, T_BASE> :: t_FitnessQuantity
            Farmer<T_GATRAITS, T_BASE> :: evaluate( t_Individual &_indiv )
            {
              if ( not _indiv.invalid() ) return t_Base::evaluate( _indiv );
              
              // fitness AND quantities of _indiv must be valid from here-on
              ++t_Base::nb_eval;
              __DODEBUGCODE( Print::out << "Evaluating Individual: " << _indiv << Print::endl; )
              unknowns.push_back( t_Unknown( &_indiv, -1) );
            }

        template<class T_GATRAITS, template< class > class T_BASE >
          void Farmer<T_GATRAITS, T_BASE> :: onWait( types::t_unsigned _bull )
          {
              Print::out << "calling next " << _bull << "." << Print::endl;
            t_Individual *indiv = next( _bull );
              Print::out << "after next " << _bull << "." << Print::endl;
            if( not indiv )
            {
              Print::out << "Sending done to bull " << _bull << "." << Print::endl;
              t_CommBase::send_command( _bull, t_CommBase::t_Commands::DONE );
              return;
            }
            Print::out << "Sending go to bull " << _bull << "." << Print::endl;
            t_CommBase::send_command( _bull, t_CommBase::t_Commands::GO );
            Print::out << "Sending individual to bull " << _bull << "." << Print::endl;
            t_CommBase::comm->send( _bull, ONWAIT_TAG( TAG ), *indiv );
            Print::out << "Activating " << _bull << "." << Print::endl;
            t_CommBase::activate( _bull );
          }

        template<class T_GATRAITS, template< class > class T_BASE >
          typename Farmer<T_GATRAITS, T_BASE> :: t_Individual*
            Farmer<T_GATRAITS, T_BASE> :: next( types::t_unsigned _bull )
            {
              namespace blamb = boost::lambda;
              typename t_Unknowns :: iterator i_unknown
                 = std::find_if( unknowns.begin(), unknowns.end(), 
                                   blamb::bind( &t_Unknown::second, blamb::_1 )
                                == blamb::constant( _bull ) );
              if ( i_unknown != unknowns.end() )
              {
                t_CommBase::comm->recv( _bull, ONWAIT_TAG( TAG+1 ), 
                                        i_unknown->first->quantities() );
                Print :: out << "received quantities: " 
                             << *(i_unknown->first) 
                             << " " << i_unknown->first->quantities() 
                             << " from " << _bull << Print::endl;
                __TRYDEBUGCODE(
                    t_Base::objective->init( *i_unknown->first );
                    i_unknown->first->set_fitness(
                      (*t_Base::objective)( i_unknown->first->const_quantities() ) );,
                    "Error while evaluating fitness.\n" )
                Print :: out << "set fitness." << Print::endl;
                unknowns.erase( i_unknown );
              }
  
              i_unknown = std::find_if( unknowns.begin(), unknowns.end(), 
                                           blamb::bind( &t_Unknown::second, blamb::_1 )
                                        == blamb::constant(-1) );
              if( i_unknown == unknowns.end() ) return NULL;
              Print :: out << "Assigning " << *i_unknown->first
                           << " to " << _bull << Print::endl;
              i_unknown->second = _bull;
              return i_unknown->first;
            }



        template<class T_GATRAITS, template< class > class T_BASE >
          void Bull<T_GATRAITS, T_BASE> :: operator()(t_Population& _parents,
                                                      t_Population& _offspring)
          {
            t_CommBase :: request( t_CommBase::t_Requests::WAITING );
            while ( t_CommBase::obey() != t_CommBase :: t_Commands :: DONE )
            {
              t_Individual individual;
              t_CommBase::comm->recv( 0, ONWAIT_TAG( TAG ), individual );
              metaeval.init( individual );
              metaeval.evaluate();
              t_CommBase::comm->send( 0, ONWAIT_TAG( TAG+1 ),
                                      individual.quantities() ); 
              t_CommBase :: request( t_CommBase::t_Requests::WAITING );
            }
            t_CommBase::command( t_CommBase::t_CowCommands::DONE );
            t_CommBase::comm->barrier();
          }

      } // namespace Evaluation
    } // namespace Graph
  } // namespace mpi
} // namespace GA
