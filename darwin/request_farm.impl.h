//
//  Version: $Id$
//
namespace GA
{
  namespace mpi
  {
    template<class T_GATRAITS>
    FarmerComm<T_GATRAITS> :: FarmerComm   ( MPI::Comm &_comm )
                                         : ::mpi::Base( _comm ),
                                           from_bulls(NULL), requests(NULL)
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

    template<class T_GATRAITS>
    FarmerComm<T_GATRAITS> :: ~FarmerComm ()
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


    template<class T_GATRAITS>
    void FarmerComm<T_GATRAITS> :: TestBulls()
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
            case REQUESTINGOBJECTIVE: derived->onObjective( *i_comp ); break;
            case REQUESTINGTABOO: derived->onTaboo( *i_comp ); break;
            case REQUESTINGHISTORYCHECK: derived->onHistory( *i_comp ); break;
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

    template< class T_GATRAITS > 
      void FarmerComm<T_GATRAITS> :: onObjective( types::t_int _bull )
      {
        
      }


    Bull :: Bull   ( MPI::Comm &_farmercomm, MPI::Comm & _cowcomm)
                 : ::mpi::Base( _farmercomm ), cowcomm(NULL),
                    to_cows( NULL ), farmer_in(t_FromFarmer::UNDEFINED),
                    farmer_out(t_ToFarmer::UNDEFINED), cows_out(NULL), ncows(0)
   {
      __ASSERT( not ::mpi::main::is_root_node() ); 
      __ASSERT( (not _cowcomm) or _cowcomm->rank() == ::mpi::ROOT_NODE ); 
   }

}


#endif
