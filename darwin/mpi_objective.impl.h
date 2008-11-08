//
//  Version: $Id$
//

namespace LaDa
{
  namespace Objective 
  {
    namespace mpi 
    {

      template<class T_GA_TRAITS >
        Farmer<T_GA_TRAITS> :: Farmer   ( MPI::Comm &_comm )
                                      : t_MPIBase( _comm ), t_Base()
        {
          try
          {
            __ASSERT( ::mpi::Base::is_root_node() ); 
          
            // Allocates Memory
            nbulls = :: mpi::Base::comm->size() -1;
            in = new types::t_unsigned[ nbulls ];
            from_bulls = new MPI::Prequest[ nbulls ];
          
            // Creates an array of persistent send requests
            MPI::Prequest *i_first = to_bulls;
            MPI::Prequest *i_end = to_bulls + nbulls;
            types::t_unsigned *i_buff = out;
            for(types::t_int i=0; i_first != i_end; ++i_first, ++i_buff, ++i )
            {
              *i_buff = t_ToBull :: UNDEFINED;
              *i_first = comm->Send_init( (const void*) i_buff, 1, mpi::INT,
                                          i, GA::mpi::Tags::Objective );
            }
          }
          catch ( MPI::exception &_e )
            __THROW_ERROR( "MPI error in Farmer constructor:\n" << _e.what() )
          catch ( std::exception &_e )
            __THROW_ERROR( "Encountered Error while constructing Farmer:\n" << _e.what() )
        }

      template< class T_GA_TRAITS >
      inline const typename Maximize<T_GA_TRAITS>::t_Fitness&
        Farmer<T_GA_TRAITS> :: evaluate( const t_Quantity &_q )
        { 
          // Sends quantities and gradients to farmer
          farmer_out = EVALUATE;
          to_farmer->start();
          to_farmer->Wait();
          // Receive results
          SendQuantity( _q, comm, 0 );
          fitness->receivefrom( comm, 0 );
          // Then broadcasts to cows
          if( not cowcomm ) return fitness();
          mpi::BroadCast bc( cowcomm );
          fitness->broadcast( _q, bc );
          bc.allocate_buffers();
          fitness->broadcast( _q, bc )
          bc.broadcast();
          fitness->broadcast( _q, bc )
          bc.reset();
          return fitness;
        }
      template< class T_GA_TRAITS >
      inline typename Maximize<T_GA_TRAITS>::t_ScalarQuantity
        Bull<T_GA_TRAITS> :: evaluate_with_gradient( const t_Quantity &_val,
                                                    t_QuantityGradients& _grad,
                                                    t_VA_Type *_i_grad)  
        { 
          types::t_unsigned n = _grad.size();
          // Sends quantities and gradients to farmer
          farmer_out = EVALUATE_WITH_GRADIENT;
          to_farmer->start();
          to_farmer->Wait();
          // Receive results
          SendQuantity( _q, comm, 0 );
          SendQuantity( _grad, comm, 0 );
          fitness->receivefrom( comm, 0 );
          comm->Recv( _i_grad, n, t_VA_Type, 0 );
          // Then broadcasts to cows
          if( not cowcomm ) return fitness();
          mpi::BroadCast bc( cowcomm );
          fitness.broadcast( bc )
          bc.serialize( _i_grad, _i_grad + n );
          bc.allocate_buffers();
          fitness.broadcast( bc )
          bc.serialize( _i_grad, _i_grad + n );
          bc.broadcast();
          fitness.broadcast( bc )
          bc.serialize( _i_grad, _i_grad + n );
          bc.reset();
          return fitness;
        }
      template< class T_GA_TRAITS >
      inline void Bull<T_GA_TRAITS> :: evaluate_gradient( const t_Quantity &_val,
                                                         t_QuantityGradients& _grad,
                                                         t_VA_Type *_i_grad)  
        { 
          types::t_unsigned n = _grad.size();
          // Sends quantities and gradients to farmer
          farmer_out = EVALUATE_GRADIENT;
          to_farmer->start();
          to_farmer->Wait();
          // Receive results
          SendQuantity( _val, comm, 0 );
          SendQuantity( _grad, comm, 0 );
          comm->Recv( _i_grad, n, t_VA_Type, 0 );
          // Then broadcasts to cows
          if( not cowcom ) return;
          mpi::BroadCast bc( cowcomm );
          bc.serialize( _i_grad, _i_grad + n );
          bc.allocate_buffers();
          bc.serialize( _i_grad, _i_grad + n );
          bc.broadcast();
          bc.serialize( _i_grad, _i_grad + n );
          bc.reset();
        }
      template< class T_GA_TRAITS >
      inline typename Maximize<T_GA_TRAITS>::t_ScalarQuantity
        Farmer<T_GA_TRAITS> :: evaluate_one_gradient( const t_Quantity &_val,
                                                      t_QuantityGradients& _grad,
                                                      types::t_unsigned _n )  
        { 
          // Sends quantities and gradients to farmer
          farmer_out = EVALUATE_ONE_GRADIENT;
          to_farmer->start();
          to_farmer->Wait();
          // Receive results
          SendQuantity( _val, comm, 0 );
          SendQuantity( _grad, comm, 0 );
          types::t_real result;
          comm->Recv( &result, 1, t_VA_Type, 0 );
          // Then broadcasts to cows
          if( not cowcom ) return;
          cowcomm->Bcast( &result, 1, t_VA_Type, 0 );
          return result;
        }


      template<class T_GA_TRAITS >
        Bull<T_GA_TRAITS> :: Bull   ( MPI::Comm &_farmercomm, MPI::Comm & _cowcomm)
                                  : t_MPIBase( _farmercomm, _cowcomm ), t_Base()
        {
           try
           {
             // Creates persistent receive/send requests with farmer
             to_farmer = new MPI::Prequest;
             *to_farmer = comm->Send_Init( (const void*) &farmer_out, 1,
                                           mpi::INT, 0, GA::mpi::Tags::Objective );
           }
           catch ( MPI::exception &_e )
             __THROW_ERROR( "MPI error in Farmer constructor:\n" << _e.what() )
           catch ( std::exception &_e )
             __THROW_ERROR( "Encountered Error while constructing Farmer:\n" << _e.what() )
        }

      template< class T_GA_TRAITS >
      inline const typename Maximize<T_GA_TRAITS>::t_Fitness&
        Bull<T_GA_TRAITS> :: evaluate( const t_Quantity &_q )
        { 
          // Sends quantities and gradients to farmer
          farmer_out = EVALUATE;
          to_farmer->start();
          to_farmer->Wait();
          // Receive results
          SendQuantity( _q, comm, 0 );
          fitness->receivefrom( comm, 0 );
          // Then broadcasts to cows
          if( not cowcom ) return fitness;
          mpi::BroadCast bc( cowcomm );
          fitness->broadcast( _q, bc );
          bc.allocate_buffers();
          fitness->broadcast( _q, bc )
          bc.broadcast();
          fitness->broadcast( _q, bc )
          bc.reset();
          return fitness;
        }
      template< class T_GA_TRAITS >
      inline typename Maximize<T_GA_TRAITS>::t_ScalarQuantity
        Bull<T_GA_TRAITS> :: evaluate_with_gradient( const t_Quantity &_val,
                                                    t_QuantityGradients& _grad,
                                                    t_VA_Type *_i_grad)  
        { 
          types::t_unsigned n = _grad.size();
          // Sends quantities and gradients to farmer
          farmer_out = EVALUATE_WITH_GRADIENT;
          to_farmer->start();
          to_farmer->Wait();
          // Receive results
          SendQuantity( _q, comm, 0 );
          SendQuantity( _grad, comm, 0 );
          fitness->receivefrom( comm, 0 );
          comm->Recv( _i_grad, n, t_VA_Type, 0 );
          // Then broadcasts to cows
          if( not cowcom ) return fitness;
          mpi::BroadCast bc( cowcomm );
          fitness.broadcast( bc )
          bc.serialize( _i_grad, _i_grad + n );
          bc.allocate_buffers();
          fitness.broadcast( bc )
          bc.serialize( _i_grad, _i_grad + n );
          bc.broadcast();
          fitness.broadcast( bc )
          bc.serialize( _i_grad, _i_grad + n );
          bc.reset();
          return fitness;
        }
      template< class T_GA_TRAITS >
      inline void Bull<T_GA_TRAITS> :: evaluate_gradient( const t_Quantity &_val,
                                                         t_QuantityGradients& _grad,
                                                         t_VA_Type *_i_grad)  
        { 
          types::t_unsigned n = _grad.size();
          // Sends quantities and gradients to farmer
          farmer_out = EVALUATE_GRADIENT;
          to_farmer->start();
          to_farmer->Wait();
          // Receive results
          SendQuantity( _val, comm, 0 );
          SendQuantity( _grad, comm, 0 );
          comm->Recv( _i_grad, n, t_VA_Type, 0 );
          // Then broadcasts to cows
          if( not cowcom ) return;
          mpi::BroadCast bc( cowcomm );
          bc.serialize( _i_grad, _i_grad + n );
          bc.allocate_buffers();
          bc.serialize( _i_grad, _i_grad + n );
          bc.broadcast();
          bc.serialize( _i_grad, _i_grad + n );
          bc.reset();
        }
      template< class T_GA_TRAITS >
      inline typename Maximize<T_GA_TRAITS>::t_ScalarQuantity
        Bull<T_GA_TRAITS> :: evaluate_one_gradient( const t_Quantity &_val,
                                                   t_QuantityGradients& _grad,
                                                   types::t_unsigned _n )  
        { 
          // Sends quantities and gradients to farmer
          farmer_out = EVALUATE_ONE_GRADIENT;
          to_farmer->start();
          to_farmer->Wait();
          // Receive results
          SendQuantity( _val, comm, 0 );
          SendQuantity( _grad, comm, 0 );
          types::t_real result;
          comm->Recv( &result, 1, t_VA_Type, 0 );
          // Then broadcasts to cows
          cowcomm->Bcast( &result, 1, t_VA_Type, 0 );
          return result;
        }
      template< class T_GA_TRAITS >
      inline const typename Maximize<T_GA_TRAITS>::t_Fitness&
        Cow<T_GA_TRAITS> :: evaluate( const t_Quantity &_q )
        { 
          mpi::BroadCast bc( *this );
          fitness->broadcast( _q, bc );
          bc.allocate_buffers();
          fitness->broadcast( _q, bc )
          bc.broadcast();
          fitness->broadcast( _q, bc )
          bc.reset();
          return fitness;
        }
      template< class T_GA_TRAITS >
      inline typename Maximize<T_GA_TRAITS>::t_ScalarQuantity
        Cow<T_GA_TRAITS> :: evaluate_with_gradient( const t_Quantity &_val,
                                                    t_QuantityGradients& _grad,
                                                    t_VA_Type *_i_grad)  
        { 
          types::t_unsigned n = _grad.size();
          mpi::BroadCast bc( *this );
          bc.serialize( _i_grad, _i_grad + n );
          fitness->broadcast( _q, bc );
          bc.allocate_buffers();
          fitness->broadcast( _q, bc )
          bc.serialize( _i_grad, _i_grad + n );
          bc.broadcast();
          fitness->broadcast( _q, bc )
          bc.serialize( _i_grad, _i_grad + n );
          bc.reset();
          return fitness;
        }
      template< class T_GA_TRAITS >
        inline void Cow<T_GA_TRAITS> :: evaluate_gradient( const t_Quantity &_val,
                                                           t_QuantityGradients& _grad,
                                                           t_VA_Type *_i_grad)  
        { 
          types::t_unsigned n = _grad.size();
          mpi::BroadCast bc( *this );
          bc.serialize( _i_grad, _i_grad + n );
          bc.allocate_buffers();
          bc.serialize( _i_grad, _i_grad + n );
          bc.broadcast();
          bc.serialize( _i_grad, _i_grad + n );
          bc.reset();
        }
      template< class T_GA_TRAITS >
      inline typename Maximize<T_GA_TRAITS>::t_ScalarQuantity
        Cow<T_GA_TRAITS> :: evaluate_one_gradient( const t_Quantity &_val,
                                                   t_QuantityGradients& _grad,
                                                   types::t_unsigned _n )  
        { 
          types::t_real result; 
          comm->Bcast( &result, 1, t_VA_Type, 0 );
          return result;
        }
    }
  }
} // namespace LaDa
