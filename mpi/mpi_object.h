//
//  Version: $Id$
//
#ifndef _MPI_OBJECT_H_
#define _MPI_OBJECT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdexcept>       // std::runtime_error
#include <math.h>
#ifdef _MPI
// #undef SEEK_SET
// #undef SEEK_CUR
// #undef SEEK_END
// #include <mpi.h>
#include <mpi2c++/mpi++.h>
#endif

#include "opt/types.h"

namespace mpi
{
#define REAL MPI::DOUBLE
#define INT MPI::INT
#define CHAR MPI::CHAR
#define UNSIGNED MPI::UNSIGNED

  extern const types::t_int ROOT_NODE;

#ifdef _MPI
  class Base
  {
    protected:
      types::t_int this_rank;
      types::t_int nproc;

    public:
      Base () : this_rank(0), nproc(1) {}
      Base   ( const Base &_c) 
           : this_rank( _c.this_rank), nproc( _c.nproc) {}
      virtual ~Base() {}

      types::t_int rank() const { return this_rank; }
      types::t_int size() const { return nproc; }
      
      bool is_root_node() const { return this_rank == ROOT_NODE; }
      void barrier() const
        { MPI::COMM_WORLD.Barrier(); }
      types::t_unsigned all_sum_all( types::t_unsigned &_in) const
      {
        types::t_unsigned out;
        MPI::COMM_WORLD.Allreduce( &_in, &out, 1, UNSIGNED, MPI::SUM ); 
        _in = out;
        return out;
      }
      types::t_int all_sum_all( types::t_int &_in) const
      {
        types::t_int out;
        MPI::COMM_WORLD.Allreduce( &_in, &out, 1, INT, MPI::SUM ); 
        _in = out;
        return out;
      }
      types::t_real all_sum_all( types::t_real &_in) const
      {
        types::t_real out;
        MPI::COMM_WORLD.Allreduce( &_in, &out, 1, REAL, MPI::SUM ); 
        _in = out;
        return out;
      }
      bool all_sum_all( bool &_bool) const
      {
        types::t_int out, in = _bool ? 1 : 0;
        MPI::COMM_WORLD.Allreduce( &in, &out, 1, UNSIGNED, MPI::SUM ); 
        _bool = ( out == nproc );
        return _bool;
      }
      bool all_or_all( bool &_bool) const
      {
        types::t_int out, in = _bool ? 1 : 0;
        MPI::COMM_WORLD.Allreduce( &in, &out, 1, UNSIGNED, MPI::SUM ); 
        _bool = ( out != 0 );
        return _bool;
      }
      bool all_xor_all( bool &_bool) const
      {
        types::t_int out, in = _bool ? 1 : 0;
        MPI::COMM_WORLD.Allreduce( &in, &out, 1, UNSIGNED, MPI::SUM ); 
        _bool = ( out == 0 or out == nproc );
        return _bool;
      }
  };

  class InitDestroy : public Base
  {
    protected:
      bool finalized;
    public:
      InitDestroy () : Base (), finalized(false) {} 
      void operator()( int _argc, char **_argv )
      { 
        if ( MPI::Is_initialized() or finalized )
          return;
        MPI::Init( _argc, _argv );
        this_rank = MPI::COMM_WORLD.Get_rank();
        nproc = MPI::COMM_WORLD.Get_size();
      }
      virtual ~InitDestroy()
      { 
        if ( ( not MPI::Is_initialized() ) or finalized )
          return;
        MPI::Finalize();
        finalized = true;
      }
  };

  class CommBase : public Base
  {
    protected:
      types::t_int *int_buff, *end_int_buff, *cur_int_buff;
      types::t_char *char_buff, *end_char_buff, *cur_char_buff;
      types::t_real *real_buff, *end_real_buff, *cur_real_buff;
      types::t_unsigned buffer_size[3];


    public:
      CommBase () : Base(), int_buff(NULL), end_int_buff(NULL), cur_int_buff(NULL),
                    char_buff(NULL), end_char_buff(NULL), cur_char_buff(NULL),
                    real_buff(NULL), end_real_buff(NULL), cur_real_buff(NULL)
      { 
        buffer_size[0] = 0;
        buffer_size[1] = 0;
        buffer_size[2] = 0;
      }
      CommBase   ( const Base &_c ) 
               : Base(_c), int_buff(NULL), end_int_buff(NULL), cur_int_buff(NULL),
                 char_buff(NULL), end_char_buff(NULL), cur_char_buff(NULL),
                 real_buff(NULL), end_real_buff(NULL), cur_real_buff(NULL)
      { 
        buffer_size[0] = 0;
        buffer_size[1] = 0;
        buffer_size[2] = 0;
      }
      CommBase   ( const CommBase &_c) 
               : Base( _c ),
                 int_buff(_c.int_buff), end_int_buff(_c.end_int_buff), cur_int_buff(_c.cur_int_buff),
                 char_buff(_c.char_buff), end_char_buff(_c.end_char_buff), cur_char_buff(_c.cur_char_buff),
                 real_buff(_c.real_buff), end_real_buff(_c.end_real_buff), cur_real_buff(_c.cur_real_buff)
      {
        buffer_size[0] = _c.buffer_size[0];
        buffer_size[1] = _c.buffer_size[1];
        buffer_size[2] = _c.buffer_size[2];
      }
      virtual ~CommBase() {}

      void destroy_buffers()
      {
        if ( int_buff ) delete[] int_buff;
        if ( char_buff ) delete[] char_buff;
        if ( real_buff ) delete[] real_buff;
        int_buff = NULL;  end_int_buff = NULL;  cur_int_buff = NULL;  buffer_size[0] = 0;  
        real_buff = NULL; end_real_buff = NULL; cur_real_buff = NULL; buffer_size[1] = 0; 
        char_buff = NULL; end_char_buff = NULL; cur_char_buff = NULL; buffer_size[2] = 0; 
      }
      virtual bool allocate_buffers( types::t_unsigned _root = ROOT_NODE ) = 0;
      virtual bool operator()( types::t_unsigned _root = ROOT_NODE ) = 0;
  };

  class BroadCast : public CommBase 
  {
    template< class T_TYPE > friend BroadCast& operator<< ( BroadCast& _this, T_TYPE &_type );
    public:
      enum t_stages { GETTING_SIZE, COPYING_TO_HERE, COPYING_FROM_HERE };
    protected:
      enum t_operation { SIZEUP, ALLOCATE, TOHERE, BROADCAST, FROMHERE, CLEAR, NEXTSTAGE };
    public:
      const static t_operation clear;
      const static t_operation nextstage;
      const static t_operation sizeup;
      const static t_operation allocate;
      const static t_operation tohere;
      const static t_operation broadcast;
      const static t_operation fromhere;

    protected:
      t_stages stage;
      bool keep_going;

    public:
      BroadCast() : CommBase(), stage(GETTING_SIZE), keep_going(true) {}
      BroadCast( const BroadCast &_b ) : CommBase(_b), stage(_b.stage) {}
      BroadCast( const Base &_b ) : CommBase(_b), stage(GETTING_SIZE) {}
      virtual ~BroadCast() { destroy_buffers(); };

      bool allocate_buffers( types::t_unsigned _root = ROOT_NODE );

      template< class T_OBJECT > bool serialize( T_OBJECT& _object );
      template< class T_ITERATOR > bool serialize( T_ITERATOR _first, T_ITERATOR _last );
      bool operator()( types::t_unsigned _root = ROOT_NODE );
      void reset()
        { destroy_buffers(); stage = GETTING_SIZE; keep_going = true; }

      t_stages get_stage() const { return stage; } 

      template< class T_CONTAINER >
      bool serialize_container ( T_CONTAINER &_cont )
      {
        types::t_int n = _cont.size();
        if ( not serialize( n ) ) return false;
        if( stage == COPYING_FROM_HERE )
          _cont.resize(n);
        typename T_CONTAINER :: iterator i_ob = _cont.begin();
        typename T_CONTAINER :: iterator i_ob_end = _cont.end();
        for(; i_ob != i_ob_end; ++i_ob )
          if( not serialize( *i_ob ) ) return false;
      
        return true;
      }
    protected:
      bool special_op( t_operation _op );
      template<class T_TYPE> bool operator_( T_TYPE &_type );
  };

  template<class T_TYPE> bool BroadCast :: operator_( T_TYPE &_type )
  {
    if ( not keep_going ) return true;
    try
    {
      if (not serialize( _type ) )
        throw std::runtime_error( "Error while BroadCasting \n" );
    }
    catch( std::exception &e )
    {
      std::cerr << "Caught error while running lada" << std::endl
                << e.what();
      keep_going = false;
      return false;
    }
    return true;
  }
  template<> inline bool BroadCast::operator_<const BroadCast::t_operation>( const t_operation &_op )
  {
    if ( not keep_going ) return true;
    return  special_op( _op );
  }

  template< class T_TYPE > inline BroadCast& operator<< ( BroadCast& _this, T_TYPE &_type )
    { _this.operator_( _type ); return _this; }

  class AllGather : public BroadCast
  {
    protected: 
      types::t_int *all_sizes;
    public:
      AllGather() : BroadCast(), all_sizes(NULL) {}
      AllGather( const AllGather &_b ) : BroadCast(_b), all_sizes(NULL) {}
      AllGather( const Base &_b ) : BroadCast(_b), all_sizes(NULL) {}
      virtual ~AllGather() 
      {
        if(all_sizes) delete[] all_sizes;
        all_sizes = NULL;
        destroy_buffers();
      };

      bool allocate_buffers( types::t_unsigned _root = ROOT_NODE );
      bool operator()( types::t_unsigned _root = ROOT_NODE );
  };


  template<class T_CONTAINER>
  inline typename T_CONTAINER::iterator begin( T_CONTAINER &_cont )
  {
    extern InitDestroy main;
    types::t_real n =   (types::t_real) main.rank() 
                      * (types::t_real) _cont.size()
                      / (types::t_real) main.size();
    return _cont.begin() + ( types::t_unsigned ) floor(n);
  }
  template<class T_CONTAINER>
  inline typename T_CONTAINER::iterator end( T_CONTAINER &_cont )
  {
    extern InitDestroy main;
    types::t_real n =   (types::t_real) ( main.rank() + 1 )
                      * (types::t_real) _cont.size()
                      / (types::t_real) main.size();
    return _cont.begin() + ( types::t_unsigned ) floor(n);
  }
  template<class T_CONTAINER>
  inline typename T_CONTAINER::const_iterator begin( const T_CONTAINER &_cont )
  {
    extern InitDestroy main;
    types::t_real n =   (types::t_real) main.rank() 
                      * (types::t_real) _cont.size()
                      / (types::t_real) main.size();
    return _cont.begin() + ( types::t_unsigned ) floor(n);
  }
  template<class T_CONTAINER>
  inline typename T_CONTAINER::const_iterator end( const T_CONTAINER &_cont )
  {
    extern InitDestroy main;
    types::t_real n =   (types::t_real) ( main.rank() + 1 )
                      * (types::t_real) _cont.size()
                      / (types::t_real) main.size();
    return _cont.begin() + ( types::t_unsigned ) floor(n);
  }
#else

  class Base
  {
     public:
       Base () {};
       Base ( const Base &_c) {};
       virtual ~Base() {}
 
       types::t_int rank() const { return 0; }
       types::t_int size() const { return 1; }
       
       void barrier() const {};
       template< class T_TYPE >
       T_TYPE all_sum_all( T_TYPE &_un ) const 
         { return _un; };
       bool all_sum_all( bool &_bool ) const 
         { return _bool; };
       bool all_or_all( bool &_bool ) const 
         { return _bool; };
       bool all_xor_all( bool &_bool ) const 
         { return _bool; };
  };

  class CommBase : virtual public Base
  {
    public:
      CommBase () {};
      CommBase ( const CommBase &_c) {};
      CommBase ( const Base &_c) {};
      virtual ~CommBase() {}

      void destroy_buffers() {};
  };
 
  class InitDestroy : virtual public Base
  {
     public:
       InitDestroy () : Base() {};
       void operator()( int, char ** ) {};
       virtual ~InitDestroy() {};
  }; 
 
  class BroadCast : virtual public CommBase 
  {
    protected:
      enum t_stages { GETTING_SIZE, COPYING_TO_HERE, COPYING_FROM_HERE };
 
    public:
      BroadCast () : CommBase () {}
      BroadCast ( const BroadCast &_b ) {}
      BroadCast ( const Base &_b )  {}
      virtual ~BroadCast() {};
 
      
      bool allocate_buffers( types::t_unsigned _root = ROOT_NODE) { return true; };
 
      template< class T_OBJECT > bool serialize( T_OBJECT _object )
        { return true; }
      template< class T_ITERATOR > bool serialize( T_ITERATOR _first, T_ITERATOR _last )
        { return true; }
      template< class T_CONTAINER > bool serialize_container ( T_CONTAINER &_cont )
        { return true; }
      bool operator()( types::t_unsigned _root = ROOT_NODE ) { return true; };
      void reset() {};
  };
  class AllGather : public BroadCast
  {
    public:
      AllGather() : BroadCast() {}
      AllGather( const AllGather &_b ) : BroadCast(_b) {}
      AllGather( const Base &_b ) : BroadCast(_b) {}
      virtual ~AllGather() {}
  };
 
  template<class T_CONTAINER>
  inline typename T_CONTAINER::iterator begin( T_CONTAINER &_cont )
    { return _cont.begin(); } 
  template<class T_CONTAINER>
  inline typename T_CONTAINER::iterator end( T_CONTAINER &_cont )
    { return _cont.end(); }

#endif


 extern InitDestroy main;


    

}
#endif
