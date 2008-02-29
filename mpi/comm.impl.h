//
//  Version: $Id$
//
namespace mpi
{
  inline CommBase::CommBase () : Base(), int_buff(NULL), end_int_buff(NULL),
                                 cur_int_buff(NULL), char_buff(NULL),
                                 end_char_buff(NULL), cur_char_buff(NULL),
                                 real_buff(NULL), end_real_buff(NULL),
                                 cur_real_buff(NULL)
  { 
    buffer_size[0] = 0;
    buffer_size[1] = 0;
    buffer_size[2] = 0;
  }
  inline CommBase::CommBase   ( MPI::Intracomm *_comm )
                            : Base( _comm ), int_buff(NULL), end_int_buff(NULL),
                              cur_int_buff(NULL), char_buff(NULL),
                              end_char_buff(NULL), cur_char_buff(NULL),
                              real_buff(NULL), end_real_buff(NULL),
                              cur_real_buff(NULL)
  { 
    buffer_size[0] = 0;
    buffer_size[1] = 0;
    buffer_size[2] = 0;
  }
  inline CommBase::CommBase   ( const Base &_c ) 
                            : Base(), int_buff(NULL), end_int_buff(NULL),
                              cur_int_buff(NULL), char_buff(NULL),
                              end_char_buff(NULL), cur_char_buff(NULL),
                              real_buff(NULL), end_real_buff(NULL),
                              cur_real_buff(NULL)
  { 
    buffer_size[0] = 0;
    buffer_size[1] = 0;
    buffer_size[2] = 0;
  }
  inline CommBase::CommBase   ( const CommBase &_c) 
                            : Base( _c ),
                              int_buff(_c.int_buff), end_int_buff(_c.end_int_buff),
                              cur_int_buff(_c.cur_int_buff),
                              char_buff(_c.char_buff), end_char_buff(_c.end_char_buff),
                              cur_char_buff(_c.cur_char_buff),
                              real_buff(_c.real_buff), end_real_buff(_c.end_real_buff),
                              cur_real_buff(_c.cur_real_buff)
  {
    buffer_size[0] = _c.buffer_size[0];
    buffer_size[1] = _c.buffer_size[1];
    buffer_size[2] = _c.buffer_size[2];
  }
  inline void CommBase::destroy_buffers()
  {
    if ( int_buff ) delete[] int_buff;
    if ( char_buff ) delete[] char_buff;
    if ( real_buff ) delete[] real_buff;
    int_buff = NULL;  end_int_buff = NULL;  cur_int_buff = NULL;  buffer_size[0] = 0;  
    real_buff = NULL; end_real_buff = NULL; cur_real_buff = NULL; buffer_size[1] = 0; 
    char_buff = NULL; end_char_buff = NULL; cur_char_buff = NULL; buffer_size[2] = 0; 
  }

  template<class T_TYPE> inline bool BroadCast :: operator_( T_TYPE &_type )
  {
    if ( not keep_going ) return true;
    __DOASSERT( not serialize( _type ),
                "Error while BroadCasting " << _type << "\n"; )
    return true;
  }

  //! Takes care of operation calls for operator<<( BroadCast&, T_TYPE&)
  template<>
    inline bool BroadCast :: 
      operator_<const BroadCast::t_operation>( const t_operation &_op )
      {
        if ( not keep_going ) return true;
        if      ( _op == SIZEUP or _op == CLEAR )   reset();
        else if ( _op == ALLOCATE ) allocate_buffers(); 
        else if ( _op == BROADCAST ) operator()(); 
        else if ( _op == NEXTSTAGE )
        {
          switch( stage )
          {
            case GETTING_SIZE: allocate_buffers(); break;
            case COPYING_TO_HERE: operator()(); break;
            case COPYING_FROM_HERE: reset(); break;
            default: break;
          }
        }
        return true; 
      }

  inline AllGather :: ~AllGather() 
  {
    if(all_sizes) delete[] all_sizes;
    all_sizes = NULL;
    destroy_buffers();
  };

  template<class T_TYPE> inline bool AllGather :: operator_( T_TYPE &_type )
  {
    if ( not keep_going ) return true;
    __DOASSERT( not serialize( _type ),
                "Error while BroadCasting " << _type << "\n"; )
    return true;
  }

  //! Takes care of operation calls for operator<<( AllGather&, T_TYPE&)
  template<> inline bool AllGather::operator_<const BroadCast::t_operation>
                                             ( const BroadCast::t_operation &_op )
  {
    if ( not keep_going ) return true;
    if      ( _op == SIZEUP or _op == CLEAR )   reset();
    else if ( _op == ALLOCATE ) allocate_buffers(); 
    else if ( _op == BROADCAST ) operator()(); 
    else if ( _op == NEXTSTAGE )
    {
      switch( stage )
      {
        case GETTING_SIZE: allocate_buffers(); break;
        case COPYING_TO_HERE: operator()(); break;
        case COPYING_FROM_HERE: reset(); break;
        default: break;
      }
    }
    return true; 
  }


  //! \brief Return an iterator so each processor gets to screen different bits
  //!        of a container
  //! \details Each proc starts at a different point in the array, and ends
  //!          right were the next begins
  //! \sa end(T_CONTAINER&), begin(T_CONTAINER&), begin(const T_CONTAINER&),
  //!     end( const T_CONTAINER&)
  //! \param _cont a container which should be the same across all processors
  template<class T_CONTAINER>
  inline typename T_CONTAINER::iterator begin( T_CONTAINER &_cont )
  {
    extern InitDestroy main;
    types::t_real n =   (types::t_real) main.rank() 
                      * (types::t_real) _cont.size()
                      / (types::t_real) main.size();
    return _cont.begin() + ( types::t_unsigned ) floor(n);
  }
  //! \brief Returns an iterator so each processor gets to screen different
  //!        bits of a container
  //! \details Each proc starts at a different point in the array, and ends
  //!          right were the next begins
  //! \sa end(T_CONTAINER&), begin(T_CONTAINER&), begin(const T_CONTAINER&),
  //!     end( const T_CONTAINER&)
  //! \param _cont a container which should be the same across all processors
  template<class T_CONTAINER>
  inline typename T_CONTAINER::iterator end( T_CONTAINER &_cont )
  {
    extern InitDestroy main;
    types::t_real n =   (types::t_real) ( main.rank() + 1 )
                      * (types::t_real) _cont.size()
                      / (types::t_real) main.size();
    return _cont.begin() + ( types::t_unsigned ) floor(n);
  }
  //! \brief Returns an constant iterator so each processor gets to screen different
  //!        bits of a container
  //! \details Each proc starts at a different point in the array, and ends
  //!          right were the next begins
  //! \sa end(T_CONTAINER&), begin(T_CONTAINER&), begin(const T_CONTAINER&),
  //!     end( const T_CONTAINER&)
  //! \param _cont a container which should be the same across all processors
  template<class T_CONTAINER>
  inline typename T_CONTAINER::const_iterator begin( const T_CONTAINER &_cont )
  {
    extern InitDestroy main;
    types::t_real n =   (types::t_real) main.rank() 
                      * (types::t_real) _cont.size()
                      / (types::t_real) main.size();
    return _cont.begin() + ( types::t_unsigned ) floor(n);
  }
  //! \brief Returns an constant iterator so each processor gets to screen different
  //!        bits of a container
  //! \details Each proc starts at a different point in the array, and ends
  //!          right were the next begins
  //! \sa end(T_CONTAINER&), begin(T_CONTAINER&), begin(const T_CONTAINER&),
  //!     end( const T_CONTAINER&)
  //! \param _cont a container which should be the same across all processors
  template<class T_CONTAINER>
  inline typename T_CONTAINER::const_iterator end( const T_CONTAINER &_cont )
  {
    extern InitDestroy main;
    types::t_real n =   (types::t_real) ( main.rank() + 1 )
                      * (types::t_real) _cont.size()
                      / (types::t_real) main.size();
    return _cont.begin() + ( types::t_unsigned ) floor(n);
  }

  //! \cond 
  template< class T_TYPE > inline AllGather& operator<< ( AllGather& _this, T_TYPE &_type )
    { _this.operator_( _type ); return _this; }
  template< class T_TYPE > inline BroadCast& operator<< ( BroadCast& _this, T_TYPE &_type )
    { _this.operator_( _type ); return _this; }
  //! \endcond

}

