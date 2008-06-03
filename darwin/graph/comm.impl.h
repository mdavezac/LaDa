//
//  Version: $Id$
//


namespace GA
{
  namespace mpi
  {
    namespace Graph
    {
      namespace Comm
      {

        template<class T_GATRAITS, class T_DERIVED>
        Farmer<T_GATRAITS, T_DERIVED> :: Farmer   ( Topology *_topo )
                                                : nbulls(0), in(NULL),
                                                  requests(NULL), taboos(NULL),
                                                  objective(NULL), store(NULL),
                                                  history(NULL),
                                                  comm( _topo->farmer_comm() )
        {
          // Allocates Memory
          nbulls = comm->size() - 1;
          in = new types::t_int[ nbulls ];
          if( not in) return;
          requests = new MPI::Prequest[ nbulls ];
          if( not requests ) { delete in; in = NULL; return; }

          // Creates an array of persistent receive requests
          MPI::Request *i_first = requests;
          MPI::Request *i_end = requests + nbulls;
          types::t_int *i_buff = in;
          for(types::t_int i=1; i_first != i_end; ++i_first, ++i_buff, ++i )
          {
            *i_buff = BullRequests :: UNDEFINED;
            *i_first = ( (MPI::Intracomm) *comm).Recv_init( i_buff, 1, MPI::INTEGER,
                                                            i, REQUEST_TAG( TAG ) );
          }
        }

        template<class T_GATRAITS, class T_DERIVED>
        Farmer<T_GATRAITS, T_DERIVED> :: ~Farmer ()
        {
          if( requests )
          {
            MPI::Prequest *i_first = requests;
            MPI::Prequest *i_end = requests + nbulls;
            for(; i_first != i_end; ++i_first ) i_first->Free();
            delete[] requests;
          }
          if( in ) delete[] in;
          requests = NULL;
          in = NULL;

          nbulls = 0;
        }

        template<class T_GATRAITS, class T_DERIVED>
        inline void Farmer<T_GATRAITS, T_DERIVED>
          :: send_command( types::t_unsigned _bull,
                           typename t_Commands :: Commands _c )
        {
          types::t_unsigned buff = _c;
          Print::out << " commanding " << _bull << " to " << buff << Print::endl;
          comm->send(_bull, COMMAND_TAG( TAG ), buff  );
        }

        template<class T_GATRAITS, class T_DERIVED>
        void Farmer<T_GATRAITS, T_DERIVED> :: test_bulls()
        {
          try
          {
            types::t_int completed[nbulls];
            types::t_int ncomp = MPI::Prequest::Waitsome( nbulls, requests,
                                                          completed );
            types::t_int *i_comp = completed;
            types::t_int *i_comp_end = completed + ncomp;
            t_Derived *derived = static_cast<t_Derived*>( this );
            for(; i_comp != i_comp_end; ++i_comp )
            {
              __ASSERT( not (*i_comp >= 0 and *i_comp < nbulls), 
                        "Process index out of range: " << *i_comp << ".\n" );

              t_Actives :: iterator i_found;
              i_found = std::find_if( actives.begin(), actives.end(),
                                      boost::lambda::_1 
                                        == boost::lambda::constant( *i_comp + 1 ) );
              __ASSERT( i_found == actives.end(),
                        "Active bull should not be active." )
              actives.erase( i_found );
#ifdef _DEBUG
              Print::out << "Bull " <<  *i_comp + 1 << " is ";
              switch( in[*i_comp] )
              {
                case t_Requests::WAITING:       Print::out << "waiting"; break;
                case t_Requests::OBJECTIVE:     Print::out << "requesting an objective"; break;
                case t_Requests::GRADIENT:      Print::out << "requesting a gradient"; break;
                case t_Requests::WITH_GRADIENT: Print::out << "requesting with gradient"; break; 
                case t_Requests::ONE_GRADIENT:  Print::out << "requesting one gradient"; break;
                case t_Requests::TABOOCHECK:    Print::out << "requesting a taboo check"; break;
                case t_Requests::HISTORYCHECK:  Print::out << "requesting a history check"; break;
                case t_Requests::STORE:         Print::out << "requesting storage"; break;
                case t_Requests::UNDEFINED:
                  __THROW_ERROR( "This request should not have been sent.\n" ) 
                  break;
              }
              Print :: out << "." << Print::endl;
#endif
              switch( in[*i_comp] )
              {
                case t_Requests::WAITING: derived->onWait( *i_comp + 1 ); break;
                case t_Requests::OBJECTIVE: derived->onObjective( *i_comp + 1 ); break;
                case t_Requests::GRADIENT: derived->onGradient( *i_comp + 1 ); break;
                case t_Requests::WITH_GRADIENT:
                  derived->onWithGradient( *i_comp + 1 ); break;
                case t_Requests::ONE_GRADIENT:
                  derived->onOneGradient( *i_comp + 1 ); break;
                case t_Requests::TABOOCHECK: derived->onTaboo( *i_comp + 1 ); break;
                case t_Requests::HISTORYCHECK: derived->onHistory( *i_comp + 1 ); break;
                case t_Requests::STORE: derived->onStore( *i_comp + 1 ); break;
                case t_Requests::UNDEFINED:
                  __THROW_ERROR( "This request should not have been sent.\n" ) 
                  break;
              }
            }
          }
          catch( std::exception &_e ) 
          {
            __THROW_ERROR(    "Encountered error while testing bulls: "
                           << _e.what() << ".\n" )
          }
          catch( MPI::Exception &_e ) 
          {
            __THROW_ERROR( "MPI error encountered.\n`" )
          }
        }

//       template<class T_GATRAITS, class T_DERIVED>
//       void Farmer<T_GATRAITS, T_DERIVED> :: wait_bulls()
//       {
//         try
//         {
//           MPI::Prequest::Waitall( nbulls, requests );
//           types::t_int *i_comp = in;
//           types::t_int *i_comp_end = in + nbulls;
//           t_Derived *derived = static_cast<t_Derived*>( this );
//           for(types::t_int bull = 1; i_comp != i_comp_end; ++i_comp, ++bull )
//           {
//             __ASSERT( not ( *i_comp >= 0 and *i_comp < nbulls ), 
//                       "Process index out of range: " << *i_comp << ".\n" );
//             switch( (typename t_Requests::Requests) *i_comp )
//             {
//               case t_Requests::WAITING:       derived->onWait( bull );         break;
//               case t_Requests::OBJECTIVE:     derived->onObjective( bull );    break;
//               case t_Requests::GRADIENT:      derived->onGradient( bull );     break;
//               case t_Requests::WITH_GRADIENT: derived->onWithGradient( bull ); break;
//               case t_Requests::ONE_GRADIENT:  derived->onOneGradient( bull );  break;
//               case t_Requests::TABOOCHECK:    derived->onTaboo( bull );        break;
//               case t_Requests::HISTORYCHECK:  derived->onHistory( bull );      break;
//               case t_Requests::STORE:         derived->onStore( bull );        break;
//               case t_Requests::UNDEFINED:
//                 __THROW_ERROR( "This request should not have been sent.\n" ) 
//                 break;
//             }
//           }
//         }
//         catch( std::exception &_e ) 
//         {
//           __THROW_ERROR(    "Encountered error while testing bulls: "
//                          << _e.what() << ".\n" )
//         }
//         catch( MPI::Exception &_e ) 
//         {
//           __THROW_ERROR( "MPI error encountered.\n`" )
//         }
//       }
        template<class T_GATRAITS, class T_DERIVED>
        inline void Farmer<T_GATRAITS, T_DERIVED> :: onTaboo( types::t_int _bull )
        {
          __ASSERT( _bull < 1 or _bull > nbulls, "bull index out of range." )
          __ASSERT( taboos, "Taboo pointer has not been set.\n")
          t_Individual indiv;
          comm->recv( _bull, ONTABOO_TAG1( TAG ), indiv );
          bool result = (*taboos)( indiv );
          comm->send( _bull, ONTABOO_TAG2( TAG ), result );
          activate(_bull);
        }
        template<class T_GATRAITS, class T_DERIVED>
        inline void Farmer<T_GATRAITS, T_DERIVED> :: onObjective( types::t_int _bull )
        {
          __ASSERT( objective, "Objective pointer not set.\n" )
          __ASSERT( _bull < 1 or _bull > nbulls, "bull index out of range." )
          t_Quantity quantities;
          t_Individual indiv;
          comm->recv( _bull, ONOBJECTIVE_TAG1( TAG ), indiv );
          objective->init( indiv );
          indiv.set_fitness( (*objective)( quantities ) );
          comm->send( _bull, ONOBJECTIVE_TAG2( TAG ), indiv.fitness() );
          activate(_bull);
        }
        template<class T_GATRAITS, class T_DERIVED>
        inline void Farmer<T_GATRAITS, T_DERIVED> :: onGradient( types::t_int _bull )
        {
          __ASSERT( objective, "Objective pointer not set.\n" )
          __ASSERT( _bull < 1 or _bull > nbulls, "bull index out of range." )
          t_Individual indiv;
          comm->recv( _bull, ONGRADIENT_TAG1( TAG ), indiv );
          t_QuantityGradients gradients;
          comm->recv( _bull, ONGRADIENT_TAG2( TAG ), gradients );
          t_VA_Type *ptrs = new t_VA_Type[ gradients.size() ];
          std::fill( ptrs, ptrs + gradients.size(), t_VA_Type(0) );
          objective->init( indiv );
          objective->evaluate_gradient( indiv.quantities(),
                                        gradients,
                                        ptrs );
          comm->send( _bull, ONGRADIENT_TAG3( TAG ), ptrs, gradients.size() );
          activate(_bull);
          delete[] ptrs;
        }
        template<class T_GATRAITS, class T_DERIVED>
        inline void Farmer<T_GATRAITS, T_DERIVED>
          :: onWithGradient( types::t_int _bull )
        {
          __ASSERT( objective, "Objective pointer not set.\n" )
          __ASSERT( _bull < 1 or _bull > nbulls, "bull index out of range." )
          t_Individual indiv;
          comm->recv( _bull, ONWITHGRADIENT_TAG1( TAG ), indiv );
          t_QuantityGradients gradients;
          comm->recv( _bull, ONWITHGRADIENT_TAG2( TAG ), gradients );
          t_VA_Type *ptrs = new t_VA_Type[ gradients.size() ];
          std::fill( ptrs, ptrs + gradients.size(), t_VA_Type(0) );
          objective->init( indiv );
          t_VA_Type result = objective->evaluate_with_gradient( indiv.quantities(),
                                                                gradients,
                                                                ptrs );
          comm->send( _bull, ONWITHGRADIENT_TAG3( TAG ), ptrs, gradients.size() );
          comm->send( _bull, ONWITHGRADIENT_TAG4( TAG ), result );
          activate(_bull);
          delete[] ptrs;
        }
        template<class T_GATRAITS, class T_DERIVED>
        inline void Farmer<T_GATRAITS, T_DERIVED>
          :: onOneGradient( types::t_int _bull )
        {
          __ASSERT( objective, "Objective pointer not set.\n" )
          __ASSERT( _bull < 1 or _bull > nbulls, "bull index out of range." )
          t_Individual indiv;
          comm->recv( _bull, ONONEGRADIENT_TAG1( TAG ), indiv );
          t_QuantityGradients gradients;
          comm->recv( _bull, ONONEGRADIENT_TAG2( TAG ), gradients );
          types::t_unsigned pos;
          comm->recv( _bull, ONONEGRADIENT_TAG3( TAG ), pos );
          types::t_real ptrs[ gradients.size() ];
          std::fill( ptrs, ptrs + gradients.size(), t_VA_Type(0) );
          objective->init( indiv );
          t_VA_Type result = objective->evaluate_one_gradient( indiv.quantities(),
                                                               gradients,
                                                               pos );
          comm->send( _bull, ONONEGRADIENT_TAG4( TAG ), result );
          activate(_bull);
        }
        template<class T_GATRAITS, class T_DERIVED>
        inline void Farmer<T_GATRAITS, T_DERIVED> :: onHistory( types::t_int _bull )
        {
          // Don't expect this message if history is not set
          __ASSERT( not history, "History Pointer not set.\n")
          __ASSERT( _bull < 1 or _bull > nbulls, "bull index out of range." )
          t_Individual indiv;
          comm->recv( _bull, ONHISTORY_TAG1( TAG ), indiv );
          bool buff = history->clone( indiv );
          comm->send( _bull, ONHISTORY_TAG2( TAG ), buff );
          activate(_bull);
          if( not buff ) return;
          comm->send( _bull, ONHISTORY_TAG3( TAG ), indiv);
        }
        template<class T_GATRAITS, class T_DERIVED>
        inline void Farmer<T_GATRAITS, T_DERIVED> :: onStore( types::t_int _bull )
        {
          __ASSERT( _bull < 1 or _bull > nbulls, "bull index out of range." )
          __ASSERT( not store, "Store Pointer not set.\n")
          t_Individual indiv;
          comm->recv( _bull, ONSTORE_TAG( TAG ), indiv );
          (*store)( indiv );
          activate(_bull);
        }
      
        template<class T_GATRAITS, class T_DERIVED>
        inline void Farmer<T_GATRAITS, T_DERIVED> :: activate( types::t_unsigned _bull )
        {
          __ASSERT( _bull < 1 or _bull > nbulls,
                    "bull index out of range." )
          actives.push_back( _bull );
          requests[_bull - 1].Start(); 
        }
        
        template<class T_GATRAITS, class T_DERIVED>
        inline void Farmer<T_GATRAITS, T_DERIVED> :: start_all()
        {
          actives.resize(nbulls);
          types::t_unsigned i = 0;
          std::generate( actives.begin(), actives.end(), ++boost::lambda::var(i) );
          MPI::Prequest::Startall( nbulls, requests ); 
        } 
      

        template<class T_GATRAITS, class T_DERIVED>
        inline void Bull<T_GATRAITS, T_DERIVED> ::
          command( const typename t_CowCommands :: Commands _c )
          {
            types::t_unsigned buff = _c;
            Print::out << " commanding " << cowcomm->size() << " cows to " << buff 
                       << " from " << cowcomm->rank() << Print::endl;
            boost::mpi::broadcast( *cowcomm, buff, 0 );
          }

        template<class T_GATRAITS, class T_DERIVED>
        inline void Bull<T_GATRAITS, T_DERIVED>
          :: request( typename t_Requests :: Requests _c ) const
        {
          types::t_unsigned buff = _c;
          MPI::Request request = ( (MPI::Intracomm) *comm).Isend( &buff, 1, 
                                                                  MPI::UNSIGNED, 0,
                                                                  REQUEST_TAG( TAG ) );
          request.Wait();
        }

        template<class T_GATRAITS, class T_DERIVED>
        inline typename Bull<T_GATRAITS, T_DERIVED> :: t_Commands :: Commands
          Bull<T_GATRAITS, T_DERIVED> :: obey()
          {
            types::t_unsigned buff;
            comm->recv( 0, COMMAND_TAG( TAG ), buff );
            Print::out << "Obeying " << buff << Print::endl; 
#ifdef _DEBUG
            switch( buff )
            {
              case 0: Print::out << "Obey: go" << Print::endl; break;
              case 1: Print::out << "Obey: done" << Print::endl; break;
              default: Print::out << "Obey: problem" << Print::endl; break;
            }
#endif
            return (typename t_Commands :: Commands) buff;
          }

        template<class T_GATRAITS, class T_DERIVED>
        inline typename Cow<T_GATRAITS, T_DERIVED> :: t_Commands :: Commands
          Cow<T_GATRAITS, T_DERIVED> :: obey()
          {
            __ASSERT( not evaluator, "Pointer to evaluator has not been set.\n" )
            types::t_unsigned buff = t_Commands::DONE;
            Print ::out << "Cow is obeying. " << comm->size() << Print::endl;
            boost::mpi::broadcast( *comm, buff, 0 );
            Print ::out << "Cow " << comm->rank() << " will ";
            t_Derived *_this = static_cast< t_Derived* >(this);
            switch( buff )
            {
              case t_Commands :: EVALUATE:
                Print::out << "evaluate." ;  break;
              case t_Commands :: WITH_GRADIENT:
                Print::out << "evaluate with gradient." ;  break;
              case t_Commands :: GRADIENT: 
                Print::out << "evaluate gradient." ;  break;
              case t_Commands :: ONE_GRADIENT:
                Print::out << "evaluate one gradient." ;  break;
              case t_Commands :: DONE: 
                Print::out << " be done." << Print::endl; break;
            }
            Print ::out << Print ::endl;
            switch( (typename t_Commands :: Commands) buff )
            {
              case t_Commands :: EVALUATE: _this->onEvaluate(); break;
              case t_Commands :: WITH_GRADIENT:
                _this->onWithGradient(); break;
              case t_Commands :: GRADIENT: _this->onGradient(); break; 
              case t_Commands :: ONE_GRADIENT:
                _this->onOneGradient(); break;
              case t_Commands :: DONE: return t_Commands :: DONE; break;
            }
            return t_Commands :: CONTINUE;
          }

        template<class T_GATRAITS, class T_DERIVED>
        inline void Cow<T_GATRAITS, T_DERIVED> :: onEvaluate()
        {
          Print :: out << "Evaluating" << Print::endl;
          typename t_GATraits :: t_Individual individual;
          Print :: out << "receiving individual" << Print::endl;
          boost::mpi::broadcast( *comm, individual, 0 );
          Print :: out << "initializing with individual" << Print::endl;
          if( not evaluator ) Print :: out << "evaluator not set " << Print::endl;
          evaluator->init( individual );
          Print :: out << "actual evaluation" << Print::endl;
          evaluator->evaluate();
          Print :: out << "Done evaluatign" << Print::endl;
        }

        template<class T_GATRAITS, class T_DERIVED>
        inline void Cow<T_GATRAITS, T_DERIVED> :: onGradient()
        {
          typename t_GATraits :: t_Individual individual;
          boost::mpi::broadcast( *comm, individual, 0 );
          gradients.resize( individual.Object().Container().size() );
          Traits::zero_out( gradients );
          evaluator->init( individual );
          evaluator->evaluate_gradient( gradients );
        }

        template<class T_GATRAITS, class T_DERIVED>
        inline void Cow<T_GATRAITS, T_DERIVED> :: onWithGradient()
        {
          typename t_GATraits :: t_Individual individual;
          boost::mpi::broadcast( *comm, individual, 0 );
          gradients.resize( individual.Object().Container().size() );
          Traits::zero_out( gradients );
          evaluator->init( individual );
          evaluator->evaluate_with_gradient( gradients );
        }

        template<class T_GATRAITS, class T_DERIVED>
        inline void Cow<T_GATRAITS, T_DERIVED> :: onOneGradient()
        {
          typename t_GATraits :: t_Individual individual;
          boost::mpi::broadcast( *comm, individual, 0 );
          types::t_unsigned pos;
          boost::mpi::broadcast( *comm, pos, 0 );
          gradients.resize( individual.Object().Container().size() );
          Traits::zero_out( gradients );
          evaluator->init( individual );
          evaluator->evaluate_one_gradient( gradients, pos );
        }

      } // namespace Comm
    } // namespace Graph

  } // namespace mpi
} // namespace GA

