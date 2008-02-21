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

        template< class T_QUANTITY >
          void send_quantities( types::t_int _bull, T_QUANTITY &_q,
                                MPI::Comm *_comm )
          {
            if( comm->rank() < 2 ) return;
            typedef Traits::Quantity< T_QUANTITY > t_QuantityTraits;
            mpi::BroadCast bc( _comm );
            t_QuantityTraits :: broadcast(_q, bc );
            bc.allocate_buffers();
            t_QuantityTraits :: broadcast(_q, bc );
            bc.send_ptp( _bull );
            t_QuantityTraits :: broadcast(_q, bc );
          }

        template< class T_QUANTITY >  
          void receive_quantities( types::t_int _bull, T_QUANTITY &_q,
                                   MPI::Comm *_comm )
          {
            if( comm->rank() < 2 ) return;
            typedef Traits::Quantity< T_QUANTITY > t_QuantityTraits;
            mpi::BroadCast bc( _comm );
            t_QuantityTraits :: broadcast(_q, bc );
            bc.allocate_buffers();
            t_QuantityTraits :: broadcast(_q, bc );
            bc.receive_ptp( _bull );
            t_QuantityTraits :: broadcast(_q, bc );
          }

        template< class T_ITERATOR >
          void send_range( types::t_int _bull, T_ITERATOR _first,
                           T_ITERATOR _last, MPI::Comm *_comm )
          {
            if( comm->rank() < 2 ) return;
            mpi::BroadCast bc( _comm );
            bc.serialize_range( _first, _last, bc );
            bc.allocate_buffers();
            bc.serialize_range( _first, _last, bc );
            bc.send_ptp( _bull );
            bc.serialize_range( _first, _last, bc );
          }

        template< class T_ITERATOR >
          void receive_range( types::t_int _bull, T_ITERATOR _first,
                              T_ITERATOR _last, MPI::Comm *_comm )
          {
            if( comm->rank() < 2 ) return;
            mpi::BroadCast bc( _comm );
            bc.serialize_range( _first, _last, bc );
            bc.allocate_buffers();
            bc.serialize_range( _first, _last, bc );
            bc.receive_ptp( _bull );
            bc.serialize_range( _first, _last, bc );
          }

        template< class T_OBJECT >  
          void bcast_template_object( types::t_int _root, T_OBJECT &_object,
                                      MPI::Comm *_comm )
          {
            if( comm->rank() < 2 ) return;
            mpi::BroadCast bc( _comm );
            _object.broadcast( bc );
            bc.allocate_buffers();
            _object.broadcast( bc ); 
            bc( _root );
            _object.broadcast( bc );
          }
        template< class T_ITERATOR >  
          void bcast_range( types::t_int _root, T_ITERATOR _first,
                            T_ITERATOR _last, MPI::Comm *_comm )
          {
            if( comm->rank() < 2 ) return;
            mpi::BroadCast bc( _comm );
            bc.serialize_range( _first, _last, bc );
            bc.allocate_buffers();
            bc.serialize_range( _first, _last, bc );
            bc( _root );
            bc.serialize_range( _first, _last, bc );
          }
        template< class T_OBJECT >  
          void bcast_object( T_OBJECT &_object, MPI::Comm *_comm )
          {
            if( comm->rank() < 2 ) return;
            mpi::BroadCast bc( _comm );
            bc << _object << mpi::BroadCast::allocate
               << _object << mpi::BroadCast::broadcast
               << _object;
          }


        template< class T_OBJECT >  
          void receive_template_object( types::t_int _bull, T_OBJECT &_object,
                                        MPI::Comm *_comm )
          {
            if( comm->rank() < 2 ) return;
            mpi::BroadCast bc( _comm );
            _object.broadcast( bc );
            bc.allocate_buffers();
            _object.broadcast( bc );
            bc.receive_ptp( _bull );
            _object.broadcast( bc );
          }

        template< class T_OBJECT >  
          void send_template_object( types::t_int _bull, T_OBJECT &_object,
                                     MPI::Comm *_comm )
          {
            if( comm->rank() < 2 ) return;
            mpi::BroadCast bc( _comm );
            _indiv.broadcast( bc );
            bc.allocate_buffers();
            _indiv.broadcast( bc );
            bc.send_ptp( _bull );
            _indiv.broadcast( bc );
          }
        template< class T_OBJECT >  
          void receive_object( types::t_int _bull,
                               T_OBJECT &_object, MPI::Comm *_comm )
          {
            if( comm->rank() < 2 ) return;
            mpi::BroadCast bc( _comm );
            bc << _object << mpi::BroadCast::allocate << _object;
            bc.receive_ptp( _bull );
            bc << _object;
          }

        template< class T_OBJECT >  
          void send_object( types::t_int _bull, T_OBJECT _object, MPI::Comm *_comm )
          {
            if( comm->rank() < 2 ) return;
            mpi::BroadCast bc( _comm );
            bc << _object << mpi::BroadCast::allocate << _object;
            bc.send_ptp( _bull );
            bc << _object;
          }

        template<class T_DERIVED>
        Farmer<T_DERIVED> :: Farmer   ( MPI::Comm *_comm )
                                    : ::mpi::Base( _comm ),
                                      from_bulls(NULL), nbulls(0),
                                      requests(NULL), taboos(NULL),
                                      objective(NULL), store(NULL),
                                      history(NULL)
        {
          __ASSERT( is_root_node() ); 

          // Allocates Memory
          nbulls = :: mpi::Base::comm->size() -1;
          in = new types::t_unsigned[ nbulls ];
          if( not in) return;
          from_bulls = new MPI::Prequest[ nbulls ];
          if( not from_bulls ) { delete in; in = NULL; return; }

          // Creates an array of persistent receive requests
          MPI::Request *i_first = from_bulls;
          MPI::Request *i_end = from_bulls + nbulls;
          types::t_int *i_buff = int;
          for(types::t_int i=0; i_first != i_end; ++i_first, ++i_buff, ++i )
          {
            *i_buff = t_FromBull :: UNDEFINED;
            *i_first = comm->Recv_init( (const void*) i_buff, 1, mpi::INT, i, TAG );
          }
        }

        template<class T_DERIVED>
        Farmer<T_DERIVED> :: ~Farmer ()
        {
          if( to_bulls )
          {
            MPI::Prequest *i_first = from_bulls;
            MPI::Prequest *i_end = from_bulls + nbulls;
            types::t_unsigned *i_buff = out;
            for(; i_first != i_end; ++i_first, ++i_buff ) i_first->Free();
            delete[] from_bulls;
          }
          if( in ) delete[] in;
          from_bulls = NULL;
          in = NULL;

          nbulls = 0;
        }

        template<class T_DERIVED>
        inline void Farmer<T_DERIVED> :: send_command( types::t_unsigned _bull,
                                                       t_Commands _c )
        {
          MPI::UNSIGNED buff = _c;
          _comm->Send( &buff, 1, MPI::UNSIGNED, _bull );
        }

        template<class T_DERIVED>
        void Farmer<T_DERIVED> :: test_bulls()
        {
          try
          {
            types::t_int *completed = types::t_int [nbulls];
            types::t_int ncomp = MPI::WaitSome( nbulls, from_bulls, completed );
            MPI::types::t_int *i_comp = completed;
            MPI::types::t_int *i_comp_end = completed + ncomp;
            t_Derived derived = static_cast<t_Derived*>( this );
            for(; i_comp != i_comp_end; ++i_comp )
            {
              __ASSERT( *i_comp > 0 and *i_comp < nbulls, 
                        "Process index out of range: " << *i_comp << ".\n" );
              switch( in[*i_comp] )
              {
                case WAITING: derived->onWait( *i_comp ); break;
                case OBJECTIVE: derived->onObjective( *i_comp ); break;
                case OBJECTIVE_GRADIENT: derived->onGradient( *i_comp ); break;
                case OBJECTIVE_WITH_GRADIENT: derived->onWithGradient( *i_comp ); break;
                case OBJECTIVE_ONE_GRADIENT: derived->onOneGradient( *i_comp ); break;
                case TABOO: derived->onTaboo( *i_comp ); break;
                case HISTORYCHECK: derived->onHistory( *i_comp ); break;
                case STIRE: derived->onStore( *i_comp ); break;
                default: __THROW_ERROR( "Unknown command" << in[*i_comp] << "." )
              }
            }

            delete[] completed;
          }
          catch( std::Exception &_e ) 
          {
            if( completed ) delete completed; completed = NULL;
            __THROW_ERROR( "Encountered error while testing bulls: " << _e.what() )
          }
          catch( MPI::Exception &_e ) 
          {
            if( completed ) delete completed; completed = NULL;
            __THROW_ERROR( "MPI error encountered: " << _e.what() )
          }
        }

        template<class T_DERIVED>
        void Farmer<T_DERIVED> :: wait_bulls()
        {
          try
          {
            types::t_int *completed = types::t_int [nbulls];
            types::t_int ncomp = MPI::WaitAll( nbulls, from_bulls, completed );
            MPI::types::t_int *i_comp = completed;
            MPI::types::t_int *i_comp_end = completed + ncomp;
            t_Derived derived = static_cast<t_Derived*>( this );
            for(; i_comp != i_comp_end; ++i_comp )
            {
              __ASSERT( *i_comp > 0 and *i_comp < nbulls, 
                        "Process index out of range: " << *i_comp << ".\n" );
              switch( in[*i_comp] )
              {
                case WAITING: derived->onWait( *i_comp ); break;
                case OBJECTIVE: derived->onObjective( *i_comp ); break;
                case OBJECTIVE_GRADIENT: derived->onGradient( *i_comp ); break;
                case OBJECTIVE_WITH_GRADIENT: derived->onWithGradient( *i_comp ); break;
                case OBJECTIVE_ONE_GRADIENT: derived->onOneGradient( *i_comp ); break;
                case TABOO: derived->onTaboo( *i_comp ); break;
                case HISTORYCHECK: derived->onHistory( *i_comp ); break;
                case STORE: derived->onStore( *i_comp ); break;
                default: __THROW_ERROR( "Unknown command" << in[*i_comp] << "." )
              }
            }

            delete[] completed;
          }
          catch( std::Exception &_e ) 
          {
            if( completed ) delete completed; completed = NULL;
            __THROW_ERROR( "Encountered error while testing bulls: " << _e.what() )
          }
          catch( MPI::Exception &_e ) 
          {
            if( completed ) delete completed; completed = NULL;
            __THROW_ERROR( "MPI error encountered: " << _e.what() )
          }
        }
        template<class T_DERIVED>
        inline void Farmer<T_DERIVED> :: onTaboo( types::t_int _bull )
        {
          __ASSERT( taboo, "Taboo pointer has not been set.\n")
          t_Individual indiv;
          t_CommBase::receive_individual( _bull, indiv );
          t_CommBase::send_object( _bull, (*taboo)( indiv ) );
          t_CommBase::activate(_bull);
        }
        template<class T_DERIVED>
        inline void Farmer<T_DERIVED> :: onObjective( types::t_int _bull )
        {
          __ASSERT( objective, "Objective pointer not set.\n" )
          typename t_Individual :: t_Quantities quantities;
          t_CommBase::receive_quantities( quantities );
          typename t_Individual :: t_Fitness fitness = (*objective)( quantities );
          t_CommBase::send_fitness( fitness );
          t_CommBase::activate(_bull);
        }
        template<class T_DERIVED>
        inline void Farmer<T_DERIVED> :: onGradient( types::t_int _bull )
        {
          __ASSERT( objective, "Objective pointer not set.\n" )
          t_Quantities quantities;
          t_CommBase::receive_quantities( quantities );
          t_QuantityGradients gradients;
          t_CommBase::receive_gradients( gradients );
          t_VA_Type *ptrs = new t_VA_Type[ gradients.size() ];
          std::fill( ptrs, ptrs + gradients.size(), t_VA_Type(0) );
          (*objective)->evaluate_gradients( quantities,
                                            gradients,
                                            ptrs );
          t_CommBase::send_gradients( _bull, ptrs, gradients.size() );
          t_CommBase::activate(_bull);
          delete[] ptrs;
        }
        template<class T_DERIVED>
        inline void Farmer<T_DERIVED> :: onWithGradient( types::t_int _bull )
        {
          __ASSERT( objective, "Objective pointer not set.\n" )
          t_Quantities quantities;
          t_CommBase::receive_quantities( quantities );
          t_QuantityGradients gradients;
          t_CommBase::receive_gradients( gradients );
          t_VA_Type *ptrs = new t_VA_Type[ gradients.size() ];
          std::fill( ptrs, ptrs + gradients.size(), t_VA_Type(0) );
          t_VA_Type result = (*objective)->evaluate_with_gradients( quantities,
                                                                    gradients,
                                                                    ptrs );
          t_CommBase::send_gradients( _bull, ptrs, gradients.size() );
          t_CommBase::send_object( _bull, result );
          t_CommBase::activate(_bull);
          delete[] ptrs;
        }
        template<class T_DERIVED>
        inline void Farmer<T_DERIVED> :: onOneGradient( types::t_int _bull )
        {
          __ASSERT( objective, "Objective pointer not set.\n" )
          t_Quantities quantities;
          t_CommBase::receive_quantities( quantities );
          t_QuantityGradients gradients;
          t_CommBase::receive_gradients( gradients );
          types::t_unsigned pos;
          t_CommBase::receive_object( pos );
          std::fill( ptrs, ptrs + gradients.size(), t_VA_Type(0) );
          t_VA_Type result = (*objective)->evaluate_one_gradients( quantities,
                                                                   gradients,
                                                                   pos );
          t_CommBase::send_object( pos );
          t_CommBase::activate(_bull);
        }
        template<class T_DERIVED>
        inline void Farmer<T_DERIVED> :: onHistory( types::t_int _bull )
        {
          // Don't expect this message if history is not set
          __ASSERT( history, "History Pointer not set.\n")
          t_Individual indiv;
          t_CommBase::receive_individual( _bull, indiv );
          MPI::BOOL buff = history->clone( indiv );
          t_CommBase::send_command( buff, _bull );
          t_CommBase::activate(_bull);
          if( not buff ) return;
          t_CommBase::send_individual( _bull, indiv );
        }
        template<class T_DERIVED>
        inline void Farmer<T_DERIVED> :: onStore( types::t_int _bull )
        {
          // Don't expect this message if history is not set
          __ASSERT( store, "History Pointer not set.\n")
          t_Individual indiv;
          t_CommBase::receive_individual( _bull, indiv );
          (*store)( indiv );
          t_CommBase::activate(_bull);
          if( not buff ) return;
          t_CommBase::send_individual( _bull, indiv );
        }
      
      

        template<class T_DERIVED>
        inline void Bull<T_DERIVED> :: command( t_Commands _c )
        {
          MPI::UNSIGNED buff = _c;
          cowcomm->Bcast( &buff, 1, MPI::UNSIGNED, 0 );
        }

        template<class T_DERIVED>
        inline t_Commands Bull<T_DERIVED> :: obey()
        {
          MPI::UNSIGNED buff;
          comm->Recv( &buff, 1, MPI::UNSIGNED, 0 );
          return buff;
        }

        template<class T_DERIVED>
        inline void Bull<T_DERIVED> :: request( t_Requests &_c )
        {
          MPI::UNSIGNED buff = _c;
          MPI::Request request = comm->Isend( &buff, 1, MPI::UNSIGNED, 0 );
          request.Wait();
        }

        template<class T_DERIVED>
        inline t_Command Cow<T_DERIVED> :: obey()
        {
          __ASSERT( evaluator, "Pointer to evaluator has not been set.\n" )
          MPI::UNSIGNED buff = t_Commands::DONE;
          comm->Bcast( &buff, 1, MPI::UNSIGNED, 0 );
          t_Derived *_this = static_cast< t_Derived* > this;
          switch( (t_Commands) buff )
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

        template<class T_DERIVED>
        inline void Cow<T_DERIVED> :: onEvaluate()
        {
          t_Individual individual;
          bcast_template_object( 0, individual, t_Base::comm );
          evaluator->init( individual );
          evaluator->evaluate();
        }

        template<class T_DERIVED>
        inline void Cow<T_DERIVED> :: onGradient()
        {
          t_Individual individual;
          bcast_template_object( 0, individual, t_Base::comm );
          gradients.resize( individual.Object().Container().size() );
          Traits::zero_out( gradients );
          evaluator->init( individual );
          evaluator->evaluate_gradients( gradients );
        }

        template<class T_DERIVED>
        inline void Cow<T_DERIVED> :: onWithGradient()
        {
          t_Individual individual;
          bcast_template_object( 0, individual, t_Base::comm );
          gradients.resize( individual.Object().Container().size() );
          Traits::zero_out( gradients );
          evaluator->init( individual );
          evaluator->evaluate_with_gradients( gradients );
        }

        template<class T_DERIVED>
        inline void Cow<T_DERIVED> :: onOneGradient()
        {
          t_Individual individual;
          bcast_template_object( 0, individual, t_Base::comm );
          MPI::UNSIGNED pos;
          bcast_object( pos, t_Base::comm );
          gradients.resize( individual.Object().Container().size() );
          Traits::zero_out( gradients );
          evaluator->init( individual );
          evaluator->evaluate_one_gradient( gradients, pos );
        }

      } // namespace Comm
    } // namespace Graph

  } // namespace mpi
} // namespace GA


#endif
