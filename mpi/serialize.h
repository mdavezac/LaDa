//
//  Version: $Id$
//
#ifndef _MPI_SERIALIZE_H_
#define _MPI_SERIALIZE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>
#include "base.h"
#ifdef ___DOTYPES
#error Cannot use ___DOTYPES
#endif
#ifdef ___DOCONTAINER 
#error Cannot use ___DOCONTAINER.
#endif
#ifdef ___TYPE__ 
#error Cannot use ___TYPE__.
#endif
#ifdef ___BUFF 
#error Cannot use ___BUFF.
#endif
#ifdef ___ENDBUFF 
#error Cannot use ___ENDBUFF.
#endif
#ifdef ___BUFFNUM
#error Cannot use ___BUFFUM.
#endif
#ifdef ___NUMTYPE
#error Cannot use ___NUMTYPE.
#endif
#ifdef ___CODE1
#error Cannot use ___CODE1.
#endif
#ifdef ___CODE2
#error Cannot use ___CODE2.
#endif


#define ___DOTYPES
#define ___NUMTYPE 0
#include "serialize.impl.h"
#define ___NUMTYPE 1 // types::t_char
#include "serialize.impl.h"
#define ___NUMTYPE 2 // types::t_real
#include "serialize.impl.h"
#define ___NUMTYPE 3 // types::t_unsigned
#include "serialize.impl.h"
#define ___NUMTYPE 4 // bool
#include "serialize.impl.h"
#undef ___DOTYPES

namespace mpi
{

  //! Object for preprocessed types
  template <class T_OBJECT> struct BroadCast::Object<const T_OBJECT>
  {
    //! Constant type
    typedef typename Modifier :: Const<T_OBJECT> :: t_constant t_constant;
    //! Non constant type
    typedef typename Modifier :: Const<T_OBJECT> :: t_nonconstant t_nonconstant;
    //! broacast instance
    BroadCast &_this;
    //! Constructor
    Object( BroadCast &__this ) : _this(__this) {}
    //! constant serializing
    bool operator()( t_constant &_ob )
     { return Object<t_nonconstant> (_this)( _ob ); }
  };

  //! Serializes range of pointers/iterator
  template <class T_OBJECT> struct BroadCast::Range
  {
    typedef typename Modifier :: Const<T_OBJECT> :: t_nonconstant t_nonconstant;
    typedef typename Modifier :: Const<T_OBJECT> :: t_constant t_constant;
    BroadCast &_this;
    Range( BroadCast &__this ) : _this(__this) {}
    //! Non constant version
    bool operator()( t_nonconstant *_first, t_nonconstant *_last )
    {
      for(; _first != _last; ++_first )
        if( not _this.serialize( *_first ) ) return false;
    
      return true;
    }
    bool operator()( t_constant *_first, t_constant *_last )
    {
      for(; _first != _last; ++_first )
        if( not _this.serialize( *_first ) ) return false;
    
      return true;
    }
  };
 
  //! \brief Contains all specialized and partially specialized classes
  //!        for serialization.
  //! \details This setup bypasses partial specialization restriction on
  //!          functions,
  template <class T_CONTAINER> struct BroadCast::Container
  {
    //! Non-constant container.
    typedef typename Modifier :: Const<T_CONTAINER> :: t_nonconstant t_nonconstant;
    //! Constant container.
    typedef typename Modifier :: Const<T_CONTAINER> :: t_constant t_constant;
    //! Reference to broadcast instance.
    BroadCast &_this;
    //! Constructor.
    Container( BroadCast &__this ) : _this(__this) {}
    //! Function for serializing
    bool operator()( t_constant &_ob );
    //! Function for serializing
    bool operator()( t_nonconstant &_ob );
  };
  template<class T_CONTAINER>
  bool BroadCast::Container<T_CONTAINER>::operator()( t_nonconstant &_cont )
  {
    types::t_unsigned n = _cont.size();
    if ( not BroadCast::Object<types::t_unsigned>(_this)( n) ) return false;
    if ( n == 0 ) return true;
    if( _this.stage == ::mpi::BroadCast::COPYING_FROM_HERE ) _cont.resize(n);
    typename T_CONTAINER :: iterator i_ob = _cont.begin();
    typename T_CONTAINER :: iterator i_ob_end = _cont.end();
    for(; i_ob != i_ob_end; ++i_ob )
      if( not _this.serialize( *i_ob ) ) return false;
  
    return true;
  }
  template<class T_CONTAINER>
  bool BroadCast::Container<T_CONTAINER>::operator()( t_constant &_cont )
  {
    if( _this.stage != COPYING_FROM_HERE ) return true;
    types::t_int n = _cont.size();
    if( not BroadCast::Object<types::t_unsigned>(_this)( n ) )
      return false;
    if ( n == 0 ) return true;
    typename T_CONTAINER :: const_iterator i_ob = _cont.begin();
    typename T_CONTAINER :: const_iterator i_ob_end = _cont.end();
    for(; i_ob != i_ob_end; ++i_ob )
      if( not _this.serialize(*i_ob ) ) return false;
    
    return true;
  };

  //! For serializing types::t_complex.
  template<> struct BroadCast::Object< types::t_complex >
  {
    BroadCast &_this;
    Object( BroadCast &__this ) : _this(__this) {}
    //! Function for serializing
    bool operator()( types::t_complex &_ob)
    {
      if( _this.stage == ::mpi::BroadCast::COPYING_TO_HERE )
      {
        if ( _this.end_real_buff - _this.cur_real_buff < 2 ) return false;
        *_this.cur_real_buff = _ob.real(); ++_this.cur_real_buff;
        *_this.cur_real_buff = _ob.imag(); ++_this.cur_real_buff;
      
        return true;
      }
      else if( _this.stage == ::mpi::BroadCast::COPYING_FROM_HERE )
      {
        if ( _this.end_real_buff - _this.cur_real_buff < 2 ) return false;
        _ob = types::t_complex( *_this.cur_real_buff, *( _this.cur_real_buff + 1 ) );
        _this.cur_real_buff += 2;
        return true;
      }
      
      _this.buffer_size[2] += 2;
      return true;
    }
    //! Function for serializing
    bool operator()( const types::t_complex &_ob )
    {
      if( _this.stage != COPYING_FROM_HERE ) return true;
      if( _this.stage == ::mpi::BroadCast::COPYING_TO_HERE )
      {
        if ( _this.end_real_buff - _this.cur_real_buff < 2 ) return false;
        *_this.cur_real_buff = _ob.real(); ++_this.cur_real_buff;
        *_this.cur_real_buff = _ob.imag(); ++_this.cur_real_buff;
      
        return true;
      }
      _this.buffer_size[2] += 2;
      return true;
    }
  };

  //! Serializes bool
  template <> struct BroadCast::Object<std::string>
  {
    //! Non constant type
    typedef Modifier :: Const<std::string> :: t_nonconstant t_nonconstant;
    //! constant type
    typedef Modifier :: Const<std::string> :: t_constant t_constant;
    //! broacast instance
    BroadCast &_this;
    //! Constructor
    Object( BroadCast &__this ) : _this(__this) {}
    //! Function for serializing
    bool operator()( t_nonconstant &_ob )
      { return Container<std::string>( _this )( _ob ); }
    //! Function for serializing
    bool operator()( t_constant &_ob )
      { return Container<std::string>( _this )( _ob ); }
  };

}
#define ___NUMTYPE 0
#define ___DOATAT 
#define ___DOCONTAINER
#include "serialize.impl.h"
#define ___NUMTYPE 2
#include "serialize.impl.h"
#undef ___DOATAT
#define ___TYPE__ types::t_char
#include "serialize.impl.h"
#define ___TYPE__ std::string
#include "serialize.impl.h"
#define ___TYPE__ types::t_complex
#include "serialize.impl.h"
#undef ___DOCONTAINER


#endif
