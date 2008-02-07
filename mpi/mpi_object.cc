//
//  Version: $Id$
//
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <stdexcept>       // std::runtime_error

#include <atat/vectmac.h>
#include <opt/debug.h>

#include "mpi_object.h"

#ifdef _MPI
namespace mpi
{
  InitDestroy main; 
  const types::t_int ROOT_NODE = 0;

  const BroadCast::t_operation BroadCast::nextstage = BroadCast::NEXTSTAGE;
  const BroadCast::t_operation BroadCast::clear     = BroadCast::CLEAR;
  const BroadCast::t_operation BroadCast::sizeup    = BroadCast::SIZEUP;
  const BroadCast::t_operation BroadCast::allocate  = BroadCast::ALLOCATE;
  const BroadCast::t_operation BroadCast::tohere    = BroadCast::TOHERE;
  const BroadCast::t_operation BroadCast::broadcast = BroadCast::BROADCAST;
  const BroadCast::t_operation BroadCast::fromhere  = BroadCast::FROMHERE;

  
  bool BroadCast :: allocate_buffers( types::t_unsigned _root )
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
    MPI::COMM_WORLD.Bcast( buffer_size, 3, UNSIGNED, _root );

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
  bool AllGather :: allocate_buffers( types::t_unsigned _root )
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
    MPI::COMM_WORLD.Allgather( buffer_size, 3, UNSIGNED, all_sizes, 3, UNSIGNED );
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
      MPI::COMM_WORLD.Bcast( int_buff, buffer_size[0], INT, _root );
    if ( char_buff and buffer_size[1] )
      MPI::COMM_WORLD.Bcast( char_buff, buffer_size[1], CHAR, _root );
    if ( real_buff and buffer_size[2] )
      MPI::COMM_WORLD.Bcast( real_buff, buffer_size[2], REAL, _root );

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
      MPI::COMM_WORLD.Allgatherv( int_buff, *(all_sizes+3*Base::this_rank), INT, int_buff,
                                  recvcounts, displs, INT );
    }

    if ( char_buff and buffer_size[1] )
    {
      for( types::t_int i=0; i < Base::nproc; ++i )
        *(recvcounts + i) = *(all_sizes + 3*i + 1);
      *displs=0;
      for( types::t_int i=1; i < Base::nproc; ++i )
        *(displs + i) = *(displs + i - 1) + *(all_sizes + 3*(i-1) + 1);
      MPI::COMM_WORLD.Allgatherv( char_buff, *(all_sizes+3*Base::this_rank+1), CHAR, char_buff,
                                  recvcounts, displs, CHAR );
    }
    if ( real_buff and buffer_size[2] )
    {
      for( types::t_int i=0; i < Base::nproc; ++i )
        *(recvcounts + i) = *(all_sizes + 3*i + 2);
      *displs=0;
      for( types::t_int i=1; i < Base::nproc; ++i )
        *(displs + i) = *(displs + i - 1) + *(all_sizes + 3*(i-1) + 2);
      MPI::COMM_WORLD.Allgatherv( real_buff, *(all_sizes+3*Base::this_rank+2), REAL, real_buff,
                                  recvcounts, displs, REAL );
    }
    
    delete[] recvcounts;
    delete[] displs;
    delete[] all_sizes; all_sizes = NULL;


    return true;
  }


  template<>
  bool BroadCast :: serialize<types::t_int>( types::t_int &_int )
  {
    if( stage != GETTING_SIZE )
    {
      if ( end_int_buff - cur_int_buff < 1) 
        return false;
      ( stage == COPYING_TO_HERE ) ? 
        *cur_int_buff = _int : _int = *cur_int_buff;
      ++cur_int_buff; 
      return true;
    }

    ++buffer_size[0]; 
    return true;
  }
  template<>
  bool BroadCast :: serialize<types::t_unsigned>( types::t_unsigned &_unsigned )
  {
    if( stage != GETTING_SIZE )
    {
      if ( end_int_buff - cur_int_buff < 1) 
        return false;
      ( stage == COPYING_TO_HERE ) ? 
        *cur_int_buff = (types::t_int) _unsigned : _unsigned = std::abs(*cur_int_buff);
      ++cur_int_buff; 
      return true;
    }

    ++buffer_size[0]; 
    return true;
  }
  template<>
  bool BroadCast :: serialize<bool>( bool &_b )
  {
    if( stage != GETTING_SIZE )
    {
      if ( end_int_buff - cur_int_buff < 1) 
        return false;
      ( stage == COPYING_TO_HERE ) ? 
        *cur_int_buff = ( _b ? 1 : 0 ) : _b =  ( *cur_int_buff == 1 ) ? 1 : 0;
      ++cur_int_buff; 
      return true;
    }

    ++buffer_size[0]; 
    return true;
  }
  template<>
  bool BroadCast :: serialize< std::vector<types::t_int> >
                             ( std::vector<types::t_int> &_cont )
  {
    if( stage == COPYING_TO_HERE )
    {
      types::t_int n = _cont.size();
      if ( end_int_buff - cur_int_buff < n+1 ) return false;
      *cur_int_buff = n; ++cur_int_buff;
      std::vector<types::t_int> :: iterator i_in = _cont.begin();
      std::vector<types::t_int> :: iterator i_end = _cont.end();
      for(; i_in != i_end; ++i_in, ++cur_real_buff )
        *cur_int_buff = *i_in;

      return true;
    }
    else if( stage == COPYING_FROM_HERE )
    {
      if ( end_int_buff - cur_int_buff < 1 ) return false;
      types::t_int n = *cur_int_buff; ++cur_int_buff;
      if ( end_int_buff - cur_int_buff < n ) return false;
      _cont.resize(n);
      std::vector<types::t_int> :: iterator i_in = _cont.begin();
      std::vector<types::t_int> :: iterator i_end = _cont.end();
      for(; i_in != i_end; ++i_in, ++cur_int_buff )
        *i_in = *cur_int_buff;

      return true;
    }

    buffer_size[0] += ( _cont.size() + 1);
    return true;
  }
  template<>
  bool BroadCast :: serialize< std::vector<types::t_unsigned> >
                             ( std::vector<types::t_unsigned> &_cont )
  {
    if( stage == COPYING_TO_HERE )
    {
      types::t_int n = _cont.size();
      if ( end_int_buff - cur_int_buff < n+1 ) return false;
      *cur_int_buff = n; ++cur_int_buff;
      std::vector<types::t_unsigned> :: iterator i_in = _cont.begin();
      std::vector<types::t_unsigned> :: iterator i_end = _cont.end();
      for(; i_in != i_end; ++i_in, ++cur_real_buff )
        *cur_int_buff = *i_in;

      return true;
    }
    else if( stage == COPYING_FROM_HERE )
    {
      if ( end_int_buff - cur_int_buff < 1 ) return false;
      types::t_int n = *cur_int_buff; ++cur_int_buff;
      if ( end_int_buff - cur_int_buff < n ) return false;
      _cont.resize(n);
      std::vector<types::t_unsigned> :: iterator i_in = _cont.begin();
      std::vector<types::t_unsigned> :: iterator i_end = _cont.end();
      for(; i_in != i_end; ++i_in, ++cur_int_buff )
        *i_in = std::abs(*cur_int_buff);

      return true;
    }

    buffer_size[0] += ( _cont.size() + 1);
    return true;
  }
  template<>
  bool BroadCast :: serialize< std::list<types::t_int> >
                             ( std::list<types::t_int> &_cont )
  {
    if( stage == COPYING_TO_HERE )
    {
      types::t_int n = _cont.size();
      if ( end_int_buff - cur_int_buff < n+1 ) return false;
      *cur_int_buff = n; ++cur_int_buff;
      std::list<types::t_int> :: iterator i_in = _cont.begin();
      std::list<types::t_int> :: iterator i_end = _cont.end();
      for(; i_in != i_end; ++i_in, ++cur_real_buff )
        *cur_int_buff = *i_in;

      return true;
    }
    else if( stage == COPYING_FROM_HERE )
    {
      if ( end_int_buff - cur_int_buff < 1 ) return false;
      types::t_int n = *cur_int_buff; ++cur_int_buff;
      if ( end_int_buff - cur_int_buff < n ) return false;
      _cont.resize(n);
      std::list<types::t_int> :: iterator i_in = _cont.begin();
      std::list<types::t_int> :: iterator i_end = _cont.end();
      for(; i_in != i_end; ++i_in, ++cur_int_buff )
        *i_in = *cur_int_buff;

      return true;
    }

    buffer_size[0] += ( _cont.size() + 1);
    return true;
  }
  template<>
  bool BroadCast :: serialize< std::list<types::t_unsigned> >
                             ( std::list<types::t_unsigned> &_cont )
  {
    if( stage == COPYING_TO_HERE )
    {
      types::t_int n = _cont.size();
      if ( end_int_buff - cur_int_buff < n+1 ) return false;
      *cur_int_buff = n; ++cur_int_buff;
      std::list<types::t_unsigned> :: iterator i_in = _cont.begin();
      std::list<types::t_unsigned> :: iterator i_end = _cont.end();
      for(; i_in != i_end; ++i_in, ++cur_real_buff )
        *cur_int_buff = *i_in;

      return true;
    }
    else if( stage == COPYING_FROM_HERE )
    {
      if ( end_int_buff - cur_int_buff < 1 ) return false;
      types::t_int n = *cur_int_buff; ++cur_int_buff;
      if ( end_int_buff - cur_int_buff < n ) return false;
      _cont.resize(n);
      std::list<types::t_unsigned> :: iterator i_in = _cont.begin();
      std::list<types::t_unsigned> :: iterator i_end = _cont.end();
      for(; i_in != i_end; ++i_in, ++cur_int_buff )
        *i_in = std::abs(*cur_int_buff);

      return true;
    }

    buffer_size[0] += ( _cont.size() + 1);
    return true;
  }

  template<>
  bool BroadCast :: serialize<std::string>( std::string &_str )
  {
    if( stage == GETTING_SIZE )
    {
      ++buffer_size[0];
      buffer_size[1] += _str.size();
      return true;
    }
    
    if ( stage == COPYING_TO_HERE )
    {
      __DOASSERT( end_int_buff - cur_int_buff < 1,
                  "Unexpected end of integer buffer\nCannot serialize string\n"; )
      std::string::iterator i_str;
      std::string::iterator i_str_end;
      *cur_int_buff = _str.size();
      ++cur_int_buff; 
      if( _str.size() == 0 ) return true;
      
      i_str = _str.begin(); i_str_end = _str.end();
      for (; i_str != i_str_end and cur_char_buff != end_char_buff; 
            ++i_str, ++cur_char_buff )
        *cur_char_buff = *i_str;
      return i_str == i_str_end;
    }

    std::string::iterator i_str;
    std::string::iterator i_str_end;
    if( end_int_buff - cur_int_buff < 1) return false;
    types::t_int size = *cur_int_buff; ++cur_int_buff; 
    _str.resize( size );
    i_str = _str.begin(); i_str_end = _str.end();
    for (; i_str != i_str_end and cur_char_buff != end_char_buff; ++i_str, ++cur_char_buff )
      *i_str =  *cur_char_buff;

    return i_str == i_str_end;
  }
  template<>
  bool BroadCast :: serialize<types::t_real>( types::t_real &_real )
  {
    if( stage != GETTING_SIZE )
    {
      if ( end_real_buff - cur_real_buff < 1 ) 
        return false;
      ( stage == COPYING_TO_HERE ) ? 
        *cur_real_buff = _real : _real = *cur_real_buff;
      ++cur_real_buff; 
      return true;
    }

    ++buffer_size[2]; 
    return true;
  }
  template<>
  bool BroadCast :: serialize<types::t_int*>( types::t_int *_first, types::t_int *_last)
  {
    if( stage != GETTING_SIZE )
    {
      if ( stage == COPYING_TO_HERE )
        for (; cur_int_buff != end_int_buff and _first != _last;
             ++cur_int_buff, ++_first )
         *cur_int_buff = *_first;
      else
        for (; cur_int_buff != end_int_buff and _first != _last;
             ++cur_int_buff, ++_first )
         *_first = *cur_int_buff;

      return _first == _last;
    }

    buffer_size[0] += (_last - _first);
    return true;
  }
  template<>
  bool BroadCast :: serialize<types::t_real*>( types::t_real *_first, types::t_real *_last)
  {
    if( stage != GETTING_SIZE )
    {
      if ( stage == COPYING_TO_HERE )
        for (; cur_real_buff != end_real_buff and _first != _last;
             ++cur_real_buff, ++_first )
         *cur_real_buff = *_first;
      else
        for (; cur_real_buff != end_real_buff and _first != _last;
             ++cur_real_buff, ++_first )
         *_first = *cur_real_buff;

      return _first == _last;
    }

    buffer_size[2] += (_last - _first);
    return true;
  }
  template<>
  bool BroadCast :: serialize< std::vector<types::t_real> >
                             ( std::vector<types::t_real> &_cont )
  {
    if( stage == COPYING_TO_HERE )
    {
      types::t_int n = _cont.size();
      if (    end_int_buff - cur_int_buff < 1 
           or end_real_buff - cur_real_buff < n ) return false;
      *cur_int_buff = n; ++cur_int_buff;
      std::vector<types::t_real> :: iterator i_in = _cont.begin();
      std::vector<types::t_real> :: iterator i_end = _cont.end();
      for(; i_in != i_end; ++i_in, ++cur_real_buff )
        *cur_real_buff = *i_in;

      return true;
    }
    else if( stage == COPYING_FROM_HERE )
    {
      if ( end_int_buff - cur_int_buff < 1 ) return false;
      types::t_int n = *cur_int_buff; ++cur_int_buff;
      if ( end_real_buff - cur_real_buff < n ) return false;
      _cont.resize(n);
      std::vector<types::t_real> :: iterator i_in = _cont.begin();
      std::vector<types::t_real> :: iterator i_end = _cont.end();
      for(; i_in != i_end; ++i_in, ++cur_real_buff )
        *i_in = *cur_real_buff;

      return true;
    }

    ++buffer_size[0];
    buffer_size[2] += _cont.size();
    return true;
  }
  template<>
  bool BroadCast :: serialize< std::list<types::t_real> >
                             ( std::list<types::t_real> &_cont )
  {
    if( stage == COPYING_TO_HERE )
    {
      types::t_int n = _cont.size();
      if (    end_int_buff - cur_int_buff < 1 
           or end_real_buff - cur_real_buff < n ) return false;
      *cur_int_buff = n; ++cur_int_buff;
      std::list<types::t_real> :: iterator i_in = _cont.begin();
      std::list<types::t_real> :: iterator i_end = _cont.end();
      for(; i_in != i_end; ++i_in, ++cur_real_buff )
        *cur_real_buff = *i_in;

      return true;
    }
    else if( stage == COPYING_FROM_HERE )
    {
      if ( end_int_buff - cur_int_buff < 1 ) return false;
      types::t_int n = *cur_int_buff; ++cur_int_buff;
      if ( end_real_buff - cur_real_buff < n ) return false;
      _cont.resize(n);
      std::list<types::t_real> :: iterator i_in = _cont.begin();
      std::list<types::t_real> :: iterator i_end = _cont.end();
      for(; i_in != i_end; ++i_in, ++cur_real_buff )
        *i_in = *cur_real_buff;

      return true;
    }

    ++buffer_size[0];
    buffer_size[2] += _cont.size();
    return true;
  }

  template<>
  bool BroadCast :: serialize< atat::rVector3d >
                             ( atat::rVector3d& _vec )
  {
    if( stage == COPYING_TO_HERE )
    {
      if ( end_real_buff - cur_real_buff < 3 ) return false;
      *cur_real_buff = _vec[0]; ++cur_real_buff;
      *cur_real_buff = _vec[1]; ++cur_real_buff;
      *cur_real_buff = _vec[2]; ++cur_real_buff;

      return true;
    }
    else if( stage == COPYING_FROM_HERE )
    {
      if ( end_real_buff - cur_real_buff < 3 ) return false;
      _vec[0] = *cur_real_buff; ++cur_real_buff;
      _vec[1] = *cur_real_buff; ++cur_real_buff;
      _vec[2] = *cur_real_buff; ++cur_real_buff;

      return true;
    }

    buffer_size[2] += 3;
    return true;
    
  }

  template<>
  bool BroadCast :: serialize< atat::rMatrix3d >
                             ( atat::rMatrix3d &_mat )
  {
    if( stage == COPYING_TO_HERE )
    {
      if ( end_real_buff - cur_real_buff < 9 ) return false;
      *cur_real_buff = _mat(0,0); ++cur_real_buff;
      *cur_real_buff = _mat(0,1); ++cur_real_buff;
      *cur_real_buff = _mat(0,2); ++cur_real_buff;
      *cur_real_buff = _mat(1,0); ++cur_real_buff;
      *cur_real_buff = _mat(1,1); ++cur_real_buff;
      *cur_real_buff = _mat(1,2); ++cur_real_buff;
      *cur_real_buff = _mat(2,0); ++cur_real_buff;
      *cur_real_buff = _mat(2,1); ++cur_real_buff;
      *cur_real_buff = _mat(2,2); ++cur_real_buff;
      return true;
    }
    if( stage == COPYING_FROM_HERE )
    {
      if ( end_real_buff - cur_real_buff < 9 ) return false;
      _mat(0,0) = *cur_real_buff; ++cur_real_buff;
      _mat(0,1) = *cur_real_buff; ++cur_real_buff;
      _mat(0,2) = *cur_real_buff; ++cur_real_buff;
      _mat(1,0) = *cur_real_buff; ++cur_real_buff;
      _mat(1,1) = *cur_real_buff; ++cur_real_buff;
      _mat(1,2) = *cur_real_buff; ++cur_real_buff;
      _mat(2,0) = *cur_real_buff; ++cur_real_buff;
      _mat(2,1) = *cur_real_buff; ++cur_real_buff;
      _mat(2,2) = *cur_real_buff; ++cur_real_buff;

      return true;
    }

    buffer_size[2] += 9;
    return true;
  }


}

#endif
