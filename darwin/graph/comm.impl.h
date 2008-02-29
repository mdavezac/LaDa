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
                                                : t_Base( *_topo->farmer_comm() ),
                                                  nbulls(0), in(NULL),
                                                  requests(NULL), taboos(NULL),
                                                  objective(NULL), store(NULL),
                                                  history(NULL)
        {
          // Allocates Memory
          nbulls = t_Base::size() - 1;
          in = new types::t_int[ nbulls ];
          if( not in) return;
          requests = new MPI::Prequest[ nbulls ];
          if( not requests ) { delete in; in = NULL; return; }

          // Creates an array of persistent receive requests
          MPI::Request *i_first = requests;
          MPI::Request *i_end = requests + nbulls;
          types::t_int *i_buff = in;
          for(types::t_int i=0; i_first != i_end; ++i_first, ++i_buff, ++i )
          {
            *i_buff = BullRequests :: UNDEFINED;
            *i_first = comm->Recv_init( i_buff, 1, MPI::INTEGER, i, TAG );
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
          t_Base::comm->Send( &buff, 1, MPI::UNSIGNED, _bull, TAG );
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
              __ASSERT( *i_comp > 0 and *i_comp < nbulls, 
                        "Process index out of range: " << *i_comp << ".\n" );
              switch( in[*i_comp] )
              {
                case t_Requests::WAITING: derived->onWait( *i_comp ); break;
                case t_Requests::OBJECTIVE: derived->onObjective( *i_comp ); break;
                case t_Requests::GRADIENT: derived->onGradient( *i_comp ); break;
                case t_Requests::WITH_GRADIENT:
                  derived->onWithGradient( *i_comp ); break;
                case t_Requests::ONE_GRADIENT:
                  derived->onOneGradient( *i_comp ); break;
                case t_Requests::TABOOCHECK: derived->onTaboo( *i_comp ); break;
                case t_Requests::HISTORYCHECK: derived->onHistory( *i_comp ); break;
                case t_Requests::STORE: derived->onStore( *i_comp ); break;
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

        template<class T_GATRAITS, class T_DERIVED>
        void Farmer<T_GATRAITS, T_DERIVED> :: wait_bulls()
        {
          try
          {
            MPI::Prequest::Waitall( nbulls, requests );
            types::t_int *i_comp = in;
            types::t_int *i_comp_end = in + nbulls;
            t_Derived *derived = static_cast<t_Derived*>( this );
            for(; i_comp != i_comp_end; ++i_comp )
            {
              __ASSERT( *i_comp > 0 and *i_comp < nbulls, 
                        "Process index out of range: " << *i_comp << ".\n" );
              switch( (typename t_Requests::Requests) *i_comp )
              {
                case t_Requests::WAITING:
                  derived->onWait( *i_comp ); break;
                case t_Requests::OBJECTIVE:
                  derived->onObjective( *i_comp ); break;
                case t_Requests::GRADIENT:
                  derived->onGradient( *i_comp ); break;
                case t_Requests::WITH_GRADIENT:
                  derived->onWithGradient( *i_comp ); break;
                case t_Requests::ONE_GRADIENT:
                  derived->onOneGradient( *i_comp ); break;
                case t_Requests::TABOOCHECK:
                  derived->onTaboo( *i_comp ); break;
                case t_Requests::HISTORYCHECK:
                  derived->onHistory( *i_comp ); break;
                case t_Requests::STORE:
                  derived->onStore( *i_comp ); break;
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
        template<class T_GATRAITS, class T_DERIVED>
        inline void Farmer<T_GATRAITS, T_DERIVED> :: onTaboo( types::t_int _bull )
        {
          __ASSERT( taboos, "Taboo pointer has not been set.\n")
          t_Individual indiv;
          receive_individual( _bull, indiv );
          send_object( _bull, (*taboos)( indiv ) );
          activate(_bull);
        }
        template<class T_GATRAITS, class T_DERIVED>
        inline void Farmer<T_GATRAITS, T_DERIVED> :: onObjective( types::t_int _bull )
        {
          __ASSERT( objective, "Objective pointer not set.\n" )
          t_Quantity quantities;
          receive_quantity( _bull, quantities );
          typename t_Individual :: t_Fitness fitness = (*objective)( quantities );
          send_fitness( _bull, fitness );
          activate(_bull);
        }
        template<class T_GATRAITS, class T_DERIVED>
        inline void Farmer<T_GATRAITS, T_DERIVED> :: onGradient( types::t_int _bull )
        {
          __ASSERT( objective, "Objective pointer not set.\n" )
          t_Quantity quantities;
          receive_quantity( _bull, quantities );
          t_QuantityGradients gradients;
          receive_gradients( _bull, gradients );
          t_VA_Type *ptrs = new t_VA_Type[ gradients.size() ];
          std::fill( ptrs, ptrs + gradients.size(), t_VA_Type(0) );
          objective->evaluate_gradient( quantities,
                                        gradients,
                                        ptrs );
          send_gradients( _bull, ptrs, gradients.size() );
          activate(_bull);
          delete[] ptrs;
        }
        template<class T_GATRAITS, class T_DERIVED>
        inline void Farmer<T_GATRAITS, T_DERIVED>
          :: onWithGradient( types::t_int _bull )
        {
          __ASSERT( objective, "Objective pointer not set.\n" )
          t_Quantity quantities;
          receive_quantity( _bull, quantities );
          t_QuantityGradients gradients;
          receive_gradients( _bull, gradients );
          t_VA_Type *ptrs = new t_VA_Type[ gradients.size() ];
          std::fill( ptrs, ptrs + gradients.size(), t_VA_Type(0) );
          t_VA_Type result = objective->evaluate_with_gradient( quantities,
                                                                gradients,
                                                                ptrs );
          send_gradients( _bull, ptrs, gradients.size() );
          send_object( _bull, result );
          activate(_bull);
          delete[] ptrs;
        }
        template<class T_GATRAITS, class T_DERIVED>
        inline void Farmer<T_GATRAITS, T_DERIVED>
          :: onOneGradient( types::t_int _bull )
        {
          __ASSERT( objective, "Objective pointer not set.\n" )
          t_Quantity quantities;
          receive_quantity( _bull, quantities );
          t_QuantityGradients gradients;
          receive_gradients( _bull, gradients );
          types::t_unsigned pos;
          receive_object( _bull, pos );
          types::t_real ptrs[ gradients.size() ];
          std::fill( ptrs, ptrs + gradients.size(), t_VA_Type(0) );
          t_VA_Type result = objective->evaluate_one_gradient( quantities,
                                                               gradients,
                                                               pos );
          send_object( _bull, pos );
          activate(_bull);
        }
        template<class T_GATRAITS, class T_DERIVED>
        inline void Farmer<T_GATRAITS, T_DERIVED> :: onHistory( types::t_int _bull )
        {
          // Don't expect this message if history is not set
          __ASSERT( history, "History Pointer not set.\n")
          t_Individual indiv;
          receive_individual( _bull, indiv );
          bool buff = history->clone( indiv );
          send_object( _bull, buff );
          activate(_bull);
          if( not buff ) return;
          send_individual( _bull, indiv );
        }
        template<class T_GATRAITS, class T_DERIVED>
        inline void Farmer<T_GATRAITS, T_DERIVED> :: onStore( types::t_int _bull )
        {
          // Don't expect this message if history is not set
          __ASSERT( store, "History Pointer not set.\n")
          t_Individual indiv;
          receive_individual( _bull, indiv );
          (*store)( indiv );
          activate(_bull);
        }
      
      

        template<class T_GATRAITS, class T_DERIVED>
        inline void Bull<T_GATRAITS, T_DERIVED> ::
          command( const typename t_CowCommands :: Commands _c )
          {
            types::t_unsigned buff = _c;
            cowcomm->Bcast( &buff, 1, MPI::UNSIGNED, 0 );
          }

        template<class T_GATRAITS, class T_DERIVED>
        inline typename Bull<T_GATRAITS, T_DERIVED> :: t_Commands :: Commands
          Bull<T_GATRAITS, T_DERIVED> :: obey()
          {
            types::t_unsigned buff;
            comm->Recv( &buff, 1, MPI::UNSIGNED, 0, TAG );
            return (typename t_Commands :: Commands) buff;
          }

        template<class T_GATRAITS, class T_DERIVED>
        inline void Bull<T_GATRAITS, T_DERIVED>
          :: request( typename t_Requests :: Requests _c ) const
        {
          types::t_unsigned buff = _c;
          MPI::Request request = comm->Isend( &buff, 1, MPI::UNSIGNED, 0, TAG );
          request.Wait();
        }

        template<class T_GATRAITS, class T_DERIVED>
        inline typename Cow<T_GATRAITS, T_DERIVED> :: t_Commands :: Commands
          Cow<T_GATRAITS, T_DERIVED> :: obey()
          {
            __ASSERT( evaluator, "Pointer to evaluator has not been set.\n" )
            types::t_unsigned buff = t_Commands::DONE;
            comm->Bcast( &buff, 1, MPI::UNSIGNED, 0 );
            t_Derived *_this = static_cast< t_Derived* >(this);
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
          typename t_GATraits :: t_Individual individual;
          ::mpi::bcast_template_object( 0, individual, t_Base::comm );
          evaluator->init( individual );
          evaluator->evaluate();
        }

        template<class T_GATRAITS, class T_DERIVED>
        inline void Cow<T_GATRAITS, T_DERIVED> :: onGradient()
        {
          typename t_GATraits :: t_Individual individual;
          ::mpi::bcast_template_object( 0, individual, t_Base::comm );
          gradients.resize( individual.Object().Container().size() );
          Traits::zero_out( gradients );
          evaluator->init( individual );
          evaluator->evaluate_gradient( gradients );
        }

        template<class T_GATRAITS, class T_DERIVED>
        inline void Cow<T_GATRAITS, T_DERIVED> :: onWithGradient()
        {
          typename t_GATraits :: t_Individual individual;
          ::mpi::bcast_template_object( 0, individual, t_Base::comm );
          gradients.resize( individual.Object().Container().size() );
          Traits::zero_out( gradients );
          evaluator->init( individual );
          evaluator->evaluate_with_gradient( gradients );
        }

        template<class T_GATRAITS, class T_DERIVED>
        inline void Cow<T_GATRAITS, T_DERIVED> :: onOneGradient()
        {
          typename t_GATraits :: t_Individual individual;
          ::mpi::bcast_template_object( 0, individual, t_Base::comm );
          types::t_unsigned pos;
          ::mpi::bcast_object( 0, pos, t_Base::comm );
          gradients.resize( individual.Object().Container().size() );
          Traits::zero_out( gradients );
          evaluator->init( individual );
          evaluator->evaluate_one_gradient( gradients, pos );
        }

      } // namespace Comm
    } // namespace Graph

  } // namespace mpi
} // namespace GA


