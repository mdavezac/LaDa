//
//  Version: $Id$
//
namespace GA
{
  namespace mpi
  {
    Farmer :: Farmer   ( MPI::Comm &_comm )
                     : ::mpi::Base( _comm ),
                       from_bulls(NULL), requests(NULL)
    {
      try
      {
        __ASSERT( ::mpi::Base::is_root_node() ); 

        // Allocates Memory
        nbulls = :: mpi::Base::comm->size() -1;
        in = new types::t_unsigned[ nbulls ];
        out = new types::t_unsigned[ nbulls ];
        to_bulls = new MPI::Prequest[ nbulls ];
        from_bulls = new MPI::Prequest[ nbulls ];

        // Creates an array of persistent send requests
        MPI::Prequest *i_first = to_bulls;
        MPI::Prequest *i_end = to_bulls + nbulls;
        types::t_unsigned *i_buff = out;
        for(types::t_int i=0; i_first != i_end; ++i_first, ++i_buff, ++i )
        {
          *i_buff = t_ToBull :: UNDEFINED;
          *i_first = comm->Send_init( (const void*) i_buff, 1, mpi::INT, i, TAG );
        }

        // Creates an array of persistent receive requests
        i_first = from_bulls;
        i_end = from_bulls + nbulls;
        i_buff = int;
        for(types::t_int i=0; i_first != i_end; ++i_first, ++i_buff, ++i )
        {
          *i_buff = t_FromBull :: UNDEFINED;
          *i_first = comm->Recv_init( (const void*) i_buff, 1, mpi::INT, i, TAG );
        }
      }
      catch ( MPI::exception &_e )
        __THROW_ERROR( "MPI error in Farmer constructor:\n" << _e.what() )
      catch ( std::exception &_e )
        __THROW_ERROR( "Encountered Error while constructing Farmer:\n" << _e.what() )
    }

    Farmer :: ~Farmer ()
    {
      try 
      {
        MPI::Prequest *i_first = to_bulls;
        MPI::Prequest *i_end = to_bulls + nbulls;
        types::t_unsigned *i_buff = out;
        if( to_bulls )
        {
          for(; i_first != i_end; ++i_first, ++i_buff ) i_first->Free();
          delete[] to_bulls;
        }
        if( out ) delete[] out;
        to_bulls = NULL;
        out = NULL;

        if( to_bulls )
        {
          i_first = from_bulls;
          i_end = from_bulls + nbulls;
          i_buff = out;
          for(; i_first != i_end; ++i_first, ++i_buff ) i_first->Free();
          delete[] from_bulls;
        }
        if( in ) delete[] in;
        from_bulls = NULL;
        in = NULL;

        nbulls = 0;
      }
      catch( MPI::Excpeption &_e )
        __THROW_ERROR( "Error encountered while freeing requests:\n" << _e.what() )
      catch( std::excpeption &_e )
        __THROW_ERROR( "Error encountered while killing farmer:\n" << _e.what() )
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
