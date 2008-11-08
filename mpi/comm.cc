//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <stdexcept>       // std::runtime_error

#include <atat/vectmac.h>
#include <opt/debug.h>

#include <print/stdout.h>

#include "mpi_object.h"

#ifdef _MPI
namespace LaDa
{
  namespace mpi
  {
    const BroadCast::t_operation BroadCast::nextstage = BroadCast::NEXTSTAGE;
    const BroadCast::t_operation BroadCast::clear     = BroadCast::CLEAR;
    const BroadCast::t_operation BroadCast::sizeup    = BroadCast::SIZEUP;
    const BroadCast::t_operation BroadCast::allocate  = BroadCast::ALLOCATE;
    const BroadCast::t_operation BroadCast::tohere    = BroadCast::TOHERE;
    const BroadCast::t_operation BroadCast::broadcast = BroadCast::BROADCAST;
    const BroadCast::t_operation BroadCast::fromhere  = BroadCast::FROMHERE;

    bool BroadCast :: allocate_buffers( types::t_unsigned _root, types::t_int _target )
    {
      if ( stage != GETTING_SIZE ) return false;

      __ASSERT(     int_buff,  "Pointer to integer buffer is already allocated.\n" )
      __ASSERT( cur_int_buff,  "Pointer to integer buffer is already allocated.\n" )
      __ASSERT( end_int_buff,  "Pointer to integer buffer is already allocated.\n" )
      __ASSERT(     char_buff,  "Pointer to character buffer is already allocated.\n" )
      __ASSERT( cur_char_buff,  "Pointer to character buffer is already allocated.\n" )
      __ASSERT( end_char_buff,  "Pointer to character buffer is already allocated.\n" )
      __ASSERT(     real_buff,  "Pointer to real-number buffer is already allocated.\n" )
      __ASSERT( cur_real_buff,  "Pointer to real-number buffer is already allocated.\n" )
      __ASSERT( end_real_buff,  "Pointer to real-number buffer is already allocated.\n" )
      stage = COPYING_TO_HERE;

      // now broadcasts sizes
      if( _target < 0 )
        comm->Bcast( buffer_size, 3, MPI::UNSIGNED, _root );
      else if( rank() == _root )
        comm->Send( buffer_size, 3, MPI::UNSIGNED, (types::t_unsigned) _target, TAG );
      else if( rank() == (types::t_unsigned) _target )
        comm->Recv( buffer_size, 3, MPI::UNSIGNED, _root, TAG );

      if ( not( buffer_size[0] or buffer_size[1] or buffer_size[2] ) )
        return false;

      // and allocates buffers
      if ( buffer_size[0]  )
      {
        int_buff = new types::t_int[ buffer_size[0] ];
        if ( not int_buff ) goto broadcast_erase;
        cur_int_buff = int_buff; end_int_buff = int_buff + buffer_size[0]; 
      }
      if ( buffer_size[1] )
      {
        char_buff = new types::t_char[ buffer_size[1] ];
        if ( not char_buff ) goto broadcast_erase;
        cur_char_buff = char_buff; end_char_buff = char_buff + buffer_size[1]; 
      }
      if ( not buffer_size[2] )
        return true;
     
      real_buff = new types::t_real[buffer_size[2]];
      if ( not real_buff ) goto broadcast_erase;
      cur_real_buff = real_buff; end_real_buff = real_buff + buffer_size[2]; 
      
      return true;

  broadcast_erase:
      if ( int_buff )  delete[] int_buff;
      int_buff = cur_int_buff = end_int_buff = NULL;
      if ( char_buff ) delete[] char_buff;
      char_buff = cur_char_buff = end_char_buff = NULL;
      if ( real_buff ) delete[] real_buff;
      real_buff = cur_real_buff = end_real_buff = NULL;
      std::cerr << "Could not allocate memory for broadcast" << std::endl;
      return false;
    }
    bool AllGather :: allocate_buffers( types::t_unsigned _root, types::t_int _target )
    {
      if ( stage != GETTING_SIZE ) return false;

      __ASSERT(     int_buff,  "Pointer to integer buffer is already allocated.\n" )
      __ASSERT( cur_int_buff,  "Pointer to integer buffer is already allocated.\n" )
      __ASSERT( end_int_buff,  "Pointer to integer buffer is already allocated.\n" )
      __ASSERT(     char_buff,  "Pointer to character buffer is already allocated.\n" )
      __ASSERT( cur_char_buff,  "Pointer to character buffer is already allocated.\n" )
      __ASSERT( end_char_buff,  "Pointer to character buffer is already allocated.\n" )
      __ASSERT(     real_buff,  "Pointer to real-number buffer is already allocated.\n" )
      __ASSERT( cur_real_buff,  "Pointer to real-number buffer is already allocated.\n" )
      __ASSERT( end_real_buff,  "Pointer to real-number buffer is already allocated.\n" )
      stage = COPYING_TO_HERE;

      // now broadcasts sizes
      all_sizes = new types::t_int[ Base::nproc * 3 ];
      comm->Allgather( buffer_size, 3, MPI::UNSIGNED, all_sizes, 3, MPI::UNSIGNED );
      buffer_size[0]=0; buffer_size[1]=0; buffer_size[2]=0;
      for( types::t_int i=0; i < Base::nproc; ++i)
      {
        buffer_size[0] += *( all_sizes + 3*i );
        buffer_size[1] += *( all_sizes + 3*i + 1 );
        buffer_size[2] += *( all_sizes + 3*i + 2 );
      }
      // and allocates buffers
      if ( buffer_size[0]  )
      {
        int_buff = new types::t_int[ buffer_size[0] ];
        if ( not int_buff ) goto gather_erase;
        cur_int_buff = int_buff; end_int_buff = int_buff + buffer_size[0]; 
      }
      if ( buffer_size[1] )
      {
        char_buff = new types::t_char[ buffer_size[1] ];
        if ( not char_buff ) goto gather_erase;
        cur_char_buff = char_buff; end_char_buff = char_buff + buffer_size[1]; 
      }
      if ( not buffer_size[2] ) return true;
     
      real_buff = new types::t_real[buffer_size[2]];
      if ( not real_buff ) goto gather_erase;
      cur_real_buff = real_buff; end_real_buff = real_buff + buffer_size[2]; 
      
      return true;

  gather_erase:
      if ( int_buff ) 
        delete[] int_buff;
      int_buff = cur_int_buff = end_int_buff = NULL;
      if ( char_buff ) 
        delete[] char_buff;
      char_buff = cur_char_buff = end_char_buff = NULL;
      if ( real_buff ) 
        delete[] real_buff;
      real_buff = cur_real_buff = end_real_buff = NULL;
      std::cerr << "Could not allocate memory for AllGatherAll" << std::endl;
      return false;
    }


    bool BroadCast :: operator()( types::t_unsigned _root )
    {
      if ( stage != COPYING_TO_HERE )
        return false;

      stage = COPYING_FROM_HERE;
      cur_int_buff = int_buff;
      cur_char_buff = char_buff;
      cur_real_buff = real_buff;
      if ( nproc == 1 )
        return true;

      if ( int_buff and buffer_size[0] )
        comm->Bcast( int_buff, buffer_size[0], MPI::INT, _root );
      if ( char_buff and buffer_size[1] )
        comm->Bcast( char_buff, buffer_size[1], MPI::CHAR, _root );
      if ( real_buff and buffer_size[2] )
        comm->Bcast( real_buff, buffer_size[2], MPI::DOUBLE, _root );

      stage = COPYING_FROM_HERE;
      cur_int_buff = int_buff;
      cur_char_buff = char_buff;
      cur_real_buff = real_buff;

      return true;
    }
    bool BroadCast :: send_ptp( types::t_unsigned _target )
    {
      if ( stage != COPYING_TO_HERE ) return false;

      stage = COPYING_FROM_HERE;
      cur_int_buff = int_buff;
      cur_char_buff = char_buff;
      cur_real_buff = real_buff;
      Print :: out << "nproc: " << nproc << Print::endl;
      if ( nproc == 1 ) return true;

      Print :: out << "buffer_size[0] " << buffer_size[0]  << Print::endl;
      Print :: out << "buffer_size[1] " << buffer_size[1]  << Print::endl;
      Print :: out << "buffer_size[2] " << buffer_size[2]  << Print::endl;
      if ( int_buff and buffer_size[0] )
        comm->Send( int_buff, buffer_size[0], MPI::INT, _target, TAG );
      if ( char_buff and buffer_size[1] )
        comm->Send( char_buff, buffer_size[1], MPI::CHAR, _target, TAG );
      if ( real_buff and buffer_size[2] )
        comm->Send( real_buff, buffer_size[2], MPI::DOUBLE, _target, TAG );

      stage = COPYING_FROM_HERE;
      cur_int_buff = int_buff;
      cur_char_buff = char_buff;
      cur_real_buff = real_buff;

      return true;
    }
    bool BroadCast :: receive_ptp( types::t_unsigned _source )
    {
      if ( stage != COPYING_TO_HERE )
        return false;

      stage = COPYING_FROM_HERE;
      cur_int_buff = int_buff;
      cur_char_buff = char_buff;
      cur_real_buff = real_buff;
      if ( nproc == 1 )
        return true;

      Print :: out << "source " << _source  << Print::endl;
      Print :: out << "buffer_size[0] " << buffer_size[0]  << Print::endl;
      if ( int_buff and buffer_size[0] )
        comm->Recv( int_buff, buffer_size[0], MPI::INT, _source, TAG );
      if ( char_buff and buffer_size[1] )
        comm->Recv( char_buff, buffer_size[1], MPI::CHAR, _source, TAG );
      if ( real_buff and buffer_size[2] )
        comm->Recv( real_buff, buffer_size[2], MPI::REAL, _source, TAG );
      Print :: out << "buffer_size[1] " << buffer_size[1]  << Print::endl;
      Print :: out << "buffer_size[2] " << buffer_size[2]  << Print::endl;

      stage = COPYING_FROM_HERE;
      cur_int_buff = int_buff;
      cur_char_buff = char_buff;
      cur_real_buff = real_buff;

      return true;
    }

    bool AllGather :: operator()( types::t_unsigned _root )
    {
      if ( stage != COPYING_TO_HERE )
        return false;

      stage = COPYING_FROM_HERE;
      cur_int_buff = int_buff;
      cur_char_buff = char_buff;
      cur_real_buff = real_buff;
      if ( nproc < 2 )
        return true;


      types::t_int *recvcounts = new types::t_int[ Base::nproc ];
      types::t_int *displs = new types::t_int[ Base::nproc ];
      if ( int_buff and buffer_size[0] )
      {
        for( types::t_int i=0; i < Base::nproc; ++i )
          *(recvcounts + i) = *(all_sizes + 3*i);
        *displs=0;
        for( types::t_int i=1; i < Base::nproc; ++i )
          *(displs + i) = *(displs + i - 1) + *(all_sizes + 3*(i-1));
        comm->Allgatherv( int_buff, *(all_sizes+3*Base::this_rank), MPI::INT,
                          int_buff, recvcounts, displs, MPI::INT );
      }

      if ( char_buff and buffer_size[1] )
      {
        for( types::t_int i=0; i < Base::nproc; ++i )
          *(recvcounts + i) = *(all_sizes + 3*i + 1);
        *displs=0;
        for( types::t_int i=1; i < Base::nproc; ++i )
          *(displs + i) = *(displs + i - 1) + *(all_sizes + 3*(i-1) + 1);
        comm->Allgatherv( char_buff, *(all_sizes+3*Base::this_rank+1), MPI::CHAR,
                          char_buff, recvcounts, displs, MPI::CHAR );
      }
      if ( real_buff and buffer_size[2] )
      {
        for( types::t_int i=0; i < Base::nproc; ++i )
          *(recvcounts + i) = *(all_sizes + 3*i + 2);
        *displs=0;
        for( types::t_int i=1; i < Base::nproc; ++i )
          *(displs + i) = *(displs + i - 1) + *(all_sizes + 3*(i-1) + 2);
        comm->Allgatherv( real_buff, *(all_sizes+3*Base::this_rank+2), MPI::DOUBLE,
                          real_buff, recvcounts, displs, MPI::DOUBLE );
      }
      
      delete[] recvcounts;
      delete[] displs;
      delete[] all_sizes; all_sizes = NULL;


      return true;
    }
  }
} // namespace LaDa
#endif
