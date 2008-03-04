//
//  Version: $Id$
//
namespace mpi
{
#include<mpi/serialize.impl.h>
  //! For serializing atat:;rVector3d
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
        *_this.___BUFF = _ob[0]; ++_this.___BUFF;
        *_this.___BUFF = _ob[1]; ++_this.___BUFF;
        *_this.___BUFF = _ob[2]; ++_this.___BUFF;
      
        return true;
      }
      else if( _this.stage == ::mpi::BroadCast::COPYING_FROM_HERE )
      {
        if ( _this.___ENDBUFF - _this.___BUFF < 3 ) return false;
        _ob[0] = *_this.___BUFF; ++_this.___BUFF;
        _ob[1] = *_this.___BUFF; ++_this.___BUFF;
        _ob[2] = *_this.___BUFF; ++_this.___BUFF;
      
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
        *_this.___BUFF = _ob[0]; ++_this.___BUFF;
        *_this.___BUFF = _ob[1]; ++_this.___BUFF;
        *_this.___BUFF = _ob[2]; ++_this.___BUFF;
      
        return true;
      }
      _this.buffer_size[___BUFFNUM] += 3;
      return true;
    }
  };
  //! For serializing atat::Matrix3d
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
        *_this.___BUFF = _ob(0,0); ++_this.___BUFF;
        *_this.___BUFF = _ob(0,1); ++_this.___BUFF;
        *_this.___BUFF = _ob(0,2); ++_this.___BUFF;
        *_this.___BUFF = _ob(1,0); ++_this.___BUFF;
        *_this.___BUFF = _ob(1,1); ++_this.___BUFF;
        *_this.___BUFF = _ob(1,2); ++_this.___BUFF;
        *_this.___BUFF = _ob(2,0); ++_this.___BUFF;
        *_this.___BUFF = _ob(2,1); ++_this.___BUFF;
        *_this.___BUFF = _ob(2,2); ++_this.___BUFF;
        return true;
      }
      if( _this.stage == ::mpi::BroadCast::COPYING_FROM_HERE )
      {
        if ( _this.___ENDBUFF - _this.___BUFF < 9 ) return false;
        _ob(0,0) = *_this.___BUFF; ++_this.___BUFF;
        _ob(0,1) = *_this.___BUFF; ++_this.___BUFF;
        _ob(0,2) = *_this.___BUFF; ++_this.___BUFF;
        _ob(1,0) = *_this.___BUFF; ++_this.___BUFF;
        _ob(1,1) = *_this.___BUFF; ++_this.___BUFF;
        _ob(1,2) = *_this.___BUFF; ++_this.___BUFF;
        _ob(2,0) = *_this.___BUFF; ++_this.___BUFF;
        _ob(2,1) = *_this.___BUFF; ++_this.___BUFF;
        _ob(2,2) = *_this.___BUFF; ++_this.___BUFF;
    
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
        *_this.___BUFF = _ob(0,0); ++_this.___BUFF;
        *_this.___BUFF = _ob(0,1); ++_this.___BUFF;
        *_this.___BUFF = _ob(0,2); ++_this.___BUFF;
        *_this.___BUFF = _ob(1,0); ++_this.___BUFF;
        *_this.___BUFF = _ob(1,1); ++_this.___BUFF;
        *_this.___BUFF = _ob(1,2); ++_this.___BUFF;
        *_this.___BUFF = _ob(2,0); ++_this.___BUFF;
        *_this.___BUFF = _ob(2,1); ++_this.___BUFF;
        *_this.___BUFF = _ob(2,2); ++_this.___BUFF;
        return true;
      }
    
      _this.buffer_size[___BUFFNUM] += 9;
      return true;
    }
  };
#endif
#ifdef ___DOATAT
  //! For serializing atat:;rVector3d
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
        *_this.___BUFF = _ob[0]; ++_this.___BUFF;
        *_this.___BUFF = _ob[1]; ++_this.___BUFF;
        *_this.___BUFF = _ob[2]; ++_this.___BUFF;
      
        return true;
      }
      else if( _this.stage == ::mpi::BroadCast::COPYING_FROM_HERE )
      {
        if ( _this.___ENDBUFF - _this.___BUFF < 3 ) return false;
        _ob[0] = *_this.___BUFF; ++_this.___BUFF;
        _ob[1] = *_this.___BUFF; ++_this.___BUFF;
        _ob[2] = *_this.___BUFF; ++_this.___BUFF;
      
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
        *_this.___BUFF = _ob[0]; ++_this.___BUFF;
        *_this.___BUFF = _ob[1]; ++_this.___BUFF;
        *_this.___BUFF = _ob[2]; ++_this.___BUFF;
      
        return true;
      }
      _this.buffer_size[___BUFFNUM] += 3;
      return true;
    }
  };
  //! For serializing atat::Matrix3d
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
        *_this.___BUFF = _ob(0,0); ++_this.___BUFF;
        *_this.___BUFF = _ob(0,1); ++_this.___BUFF;
        *_this.___BUFF = _ob(0,2); ++_this.___BUFF;
        *_this.___BUFF = _ob(1,0); ++_this.___BUFF;
        *_this.___BUFF = _ob(1,1); ++_this.___BUFF;
        *_this.___BUFF = _ob(1,2); ++_this.___BUFF;
        *_this.___BUFF = _ob(2,0); ++_this.___BUFF;
        *_this.___BUFF = _ob(2,1); ++_this.___BUFF;
        *_this.___BUFF = _ob(2,2); ++_this.___BUFF;
        return true;
      }
      if( _this.stage == ::mpi::BroadCast::COPYING_FROM_HERE )
      {
        if ( _this.___ENDBUFF - _this.___BUFF < 9 ) return false;
        _ob(0,0) = *_this.___BUFF; ++_this.___BUFF;
        _ob(0,1) = *_this.___BUFF; ++_this.___BUFF;
        _ob(0,2) = *_this.___BUFF; ++_this.___BUFF;
        _ob(1,0) = *_this.___BUFF; ++_this.___BUFF;
        _ob(1,1) = *_this.___BUFF; ++_this.___BUFF;
        _ob(1,2) = *_this.___BUFF; ++_this.___BUFF;
        _ob(2,0) = *_this.___BUFF; ++_this.___BUFF;
        _ob(2,1) = *_this.___BUFF; ++_this.___BUFF;
        _ob(2,2) = *_this.___BUFF; ++_this.___BUFF;
    
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
        *_this.___BUFF = _ob(0,0); ++_this.___BUFF;
        *_this.___BUFF = _ob(0,1); ++_this.___BUFF;
        *_this.___BUFF = _ob(0,2); ++_this.___BUFF;
        *_this.___BUFF = _ob(1,0); ++_this.___BUFF;
        *_this.___BUFF = _ob(1,1); ++_this.___BUFF;
        *_this.___BUFF = _ob(1,2); ++_this.___BUFF;
        *_this.___BUFF = _ob(2,0); ++_this.___BUFF;
        *_this.___BUFF = _ob(2,1); ++_this.___BUFF;
        *_this.___BUFF = _ob(2,2); ++_this.___BUFF;
        return true;
      }
    
      _this.buffer_size[___BUFFNUM] += 9;
      return true;
    }
  };
}
