//
//  Version: $Id$
//

#ifndef _SERIALIZE_ATAT_H_

#ifndef ___BUFFNUM

#define ___BUFFNUM 0
#include "serialize.h"
#define ___BUFFNUM 2
#include "serialize.h"
#undef ___BUFFNUM

#define _SERIALIZE_ATAT_H_

#else

#if ___BUFFNUM == 0
#define ___TYPE__  types::t_int
#define ___BUFF    cur_int_buff
#define ___ENDBUFF end_int_buff
#define ___CODE1( code ) code
#define ___CODE2( code ) code
#elif ___BUFFNUM == 1 
#define ___TYPE__  types::t_char
#define ___BUFF    cur_char_buff
#define ___ENDBUFF end_char_buff
#define ___CODE1( code ) code
#define ___CODE2( code ) code
#elif ___BUFFNUM == 2 
#define ___TYPE__  types::t_real
#define ___BUFF    cur_real_buff
#define ___ENDBUFF end_real_buff
#define ___CODE1( code ) code
#define ___CODE2( code ) code
#elif ___BUFFNUM == 3
#define ___TYPE__  types::t_unsigned
#define ___BUFF    cur_int_buff
#define ___ENDBUFF end_int_buff
#define ___CODE1( code ) (___TYPE__) code 
#define ___CODE2( code ) std::abs(code)
#elif ___BUFFNUM == 4
#define ___TYPE__  bool
#define ___BUFF    cur_int_buff
#define ___ENDBUFF end_int_buff
#define ___CODE1( code ) code ? 1: 0
#define ___CODE2( code ) code == 0
#endif

namespace mpi
{
  /** \ingroup MPI
  * \brief Serializes an  atat:;rVector3d
  */
  template<> struct BroadCast::Object< atat::Vector3d<___TYPE__> >
  {
    //! Type of the object to serialize.
    typedef atat::Vector3d<___TYPE__> t_Object;
    //! Non constant type
    typedef Modifier :: Const<t_Object> :: t_nonconstant t_nonconstant;
    //! constant type
    typedef Modifier :: Const<t_Object> :: t_constant t_constant;
    //! broacast instance
    BroadCast &_this;
    //! Constructor
    Object( BroadCast &__this ) : _this(__this) {}
    //! Function for serializing
    bool operator()( t_nonconstant &_ob )
    {
      if( _this.stage == ::mpi::BroadCast::COPYING_TO_HERE )
      {
        if ( _this.___ENDBUFF - _this.___BUFF < 3 ) return false;
        *_this.___BUFF = ___CODE1(_ob[0]); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob[1]); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob[2]); ++_this.___BUFF;
      
        return true;
      }
      else if( _this.stage == ::mpi::BroadCast::COPYING_FROM_HERE )
      {
        if ( _this.___ENDBUFF - _this.___BUFF < 3 ) return false;
        _ob[0] = ___CODE2(*_this.___BUFF); ++_this.___BUFF;
        _ob[1] = ___CODE2(*_this.___BUFF); ++_this.___BUFF;
        _ob[2] = ___CODE2(*_this.___BUFF); ++_this.___BUFF;
      
        return true;
      }
      
      _this.buffer_size[___BUFFNUM] += 3;
      return true;
    }
    //! Function for serializing
    bool operator()( t_constant &_ob )
    {
      __ASSERT( _this.stage == COPYING_FROM_HERE,
                "Attempting to change a constant object.\n" )
      if( _this.stage == ::mpi::BroadCast::COPYING_TO_HERE )
      {
        if ( _this.___ENDBUFF - _this.___BUFF < 3 ) return false;
        *_this.___BUFF = ___CODE1(_ob[0]); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob[1]); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob[2]); ++_this.___BUFF;
      
        return true;
      }
      _this.buffer_size[___BUFFNUM] += 3;
      return true;
    }
  };
  /** \ingroup MPI
  * \brief Serializes an  atat:;rVector3d
  */
  template<> struct BroadCast::Object< atat::Matrix3d< ___TYPE__ > >
  {
    //! Type of the object to serialize.
    typedef atat::Matrix3d<___TYPE__> t_Object;
    //! Non constant type
    typedef Modifier :: Const<t_Object> :: t_nonconstant t_nonconstant;
    //! constant type
    typedef Modifier :: Const<t_Object> :: t_constant t_constant;
    //! broacast instance
    BroadCast &_this;
    //! Constructor
    Object( BroadCast &__this ) : _this(__this) {}
    //! Function for serializing
    bool operator()( t_nonconstant &_ob )
    {
      if( _this.stage == ::mpi::BroadCast::COPYING_TO_HERE )
      {
        if ( _this.___ENDBUFF - _this.___BUFF < 9 ) return false;
        *_this.___BUFF = ___CODE1(_ob(0,0)); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob(0,1)); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob(0,2)); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob(1,0)); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob(1,1)); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob(1,2)); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob(2,0)); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob(2,1)); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob(2,2)); ++_this.___BUFF;
        return true;
      }
      if( _this.stage == ::mpi::BroadCast::COPYING_FROM_HERE )
      {
        if ( _this.___ENDBUFF - _this.___BUFF < 9 ) return false;
        _ob(0,0) = ___CODE2(*_this.___BUFF); ++_this.___BUFF;
        _ob(0,1) = ___CODE2(*_this.___BUFF); ++_this.___BUFF;
        _ob(0,2) = ___CODE2(*_this.___BUFF); ++_this.___BUFF;
        _ob(1,0) = ___CODE2(*_this.___BUFF); ++_this.___BUFF;
        _ob(1,1) = ___CODE2(*_this.___BUFF); ++_this.___BUFF;
        _ob(1,2) = ___CODE2(*_this.___BUFF); ++_this.___BUFF;
        _ob(2,0) = ___CODE2(*_this.___BUFF); ++_this.___BUFF;
        _ob(2,1) = ___CODE2(*_this.___BUFF); ++_this.___BUFF;
        _ob(2,2) = ___CODE2(*_this.___BUFF); ++_this.___BUFF;
    
        return true;
      }
    
      _this.buffer_size[___BUFFNUM] += 9;
      return true;
    }
    //! Function for serializing
    bool operator()( t_constant &_ob )
    {
      __ASSERT( _this.stage == COPYING_FROM_HERE,
                "Attempting to change a constant object.\n" )
      if( _this.stage == ::mpi::BroadCast::COPYING_TO_HERE )
      {
        if ( _this.___ENDBUFF - _this.___BUFF < 9 ) return false;
        *_this.___BUFF = ___CODE1(_ob(0,0)); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob(0,1)); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob(0,2)); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob(1,0)); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob(1,1)); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob(1,2)); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob(2,0)); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob(2,1)); ++_this.___BUFF;
        *_this.___BUFF = ___CODE1(_ob(2,2)); ++_this.___BUFF;
        return true;
      }
    
      _this.buffer_size[___BUFFNUM] += 9;
      return true;
    }
  };
} // namespace atat

#undef ___TYPE__
#undef ___BUFF
#undef ___ENDBUFF
#undef ___BUFFNUM
#undef ___CODE1
#undef ___CODE2

#endif
#endif
