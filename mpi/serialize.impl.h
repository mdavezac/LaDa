//
//  Version: $Id$
//
#ifndef ___TYPE__
#if ___NUMTYPE == 0
#define ___BUFFNUM 0
#define ___TYPE__  types::t_int
#define ___BUFF    cur_int_buff
#define ___ENDBUFF end_int_buff
#define ___CODE1( code ) code
#define ___CODE2( code ) code
#elif ___NUMTYPE == 1 
#define ___BUFFNUM 1
#define ___TYPE__  types::t_char
#define ___BUFF    cur_char_buff
#define ___ENDBUFF end_char_buff
#define ___CODE1( code ) code
#define ___CODE2( code ) code
#elif ___NUMTYPE == 2 
#define ___BUFFNUM 2
#define ___TYPE__  types::t_real
#define ___BUFF    cur_real_buff
#define ___ENDBUFF end_real_buff
#define ___CODE1( code ) code
#define ___CODE2( code ) code
#elif ___NUMTYPE == 3
#define ___BUFFNUM 0
#define ___TYPE__  types::t_unsigned
#define ___BUFF    cur_int_buff
#define ___ENDBUFF end_int_buff
#define ___CODE1( code ) (___TYPE__) code 
#define ___CODE2( code ) std::abs(code)
#elif ___NUMTYPE == 4
#define ___BUFFNUM 0
#define ___TYPE__  bool
#define ___BUFF    cur_int_buff
#define ___ENDBUFF end_int_buff
#define ___CODE1( code ) code ? 1: 0
#define ___CODE2( code ) code == 1
#endif
#endif

namespace mpi
{
#ifdef ___DOTYPES
  //! Object for preprocessed types
  template <> struct BroadCast::Object<___TYPE__>
  {
    //! Non constant type
    typedef Modifier :: Const<___TYPE__> :: t_nonconstant t_nonconstant;
    //! constant type
    typedef Modifier :: Const<___TYPE__> :: t_constant t_constant;
    //! broacast instance
    BroadCast &_this;
    //! Constructor
    Object( BroadCast &__this ) : _this(__this) {}
    //! non constant serializing
    bool operator()( t_nonconstant &_ob )
    { 
      if( _this.stage == ::mpi::BroadCast::COPYING_TO_HERE )
      {
        if ( _this.___ENDBUFF - _this.___BUFF < 1 ) return false;
        *_this.___BUFF = ___CODE1(_ob); ++_this.___BUFF;
      
        return true;
      }
      else if( _this.stage == ::mpi::BroadCast::COPYING_FROM_HERE )
      {
        if ( _this.___ENDBUFF - _this.___BUFF < 1 ) return false;
         _ob = ___CODE2(*_this.___BUFF); ++_this.___BUFF;
        return true;
      }
      
//     types::t_int a = ___BUFFNUM;
      _this.buffer_size[___BUFFNUM] += 1;
      return true;
    }
    //! constant serializing
    bool operator()( t_constant &_ob )
    { 
      if( _this.stage == COPYING_FROM_HERE ) return true;
      if( _this.stage == ::mpi::BroadCast::COPYING_TO_HERE )
      {
        if ( _this.___ENDBUFF - _this.___BUFF < 1 ) return false;
        *_this.___BUFF = ___CODE1(_ob); ++_this.___BUFF;
      
        return true;
      }
      _this.buffer_size[___BUFFNUM] += 1;
      return true;
    }
  };
#endif

#ifdef ___DOCONTAINER
  //! Serializes vectors
  template <> struct BroadCast::Object< std::list< ___TYPE__ > >
  {
    //! type of the object
    typedef std::vector<___TYPE__> t_Object;
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
      { return Container<t_Object>( _this )( _ob ); }
    //! Function for serializing
    bool operator()( t_constant &_ob )
      { return Container<t_Object>( _this )( _ob ); }
  };
  //! Serializes vectors
  template <> struct BroadCast::Object< std::vector< ___TYPE__ > >
  {
    //! type of the object
    typedef std::vector<___TYPE__> t_Object;
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
      { return Container<t_Object>( _this )( _ob ); }
    //! Function for serializing
    bool operator()( t_constant &_ob )
      { return Container<t_Object>( _this )( _ob ); }
  };
#endif

}
#if defined( ___OBJECTCODE )
#if defined( ___TYPE__ )
  template <> struct BroadCast::Object< ___TYPE__ >
  {
    //! Non constant type
    typedef ___TYPE__ t_Object;
    //! Non constant type
    typedef Modifier :: Const<t_Object> :: t_nonconstant t_nonconstant;
    //! constant type
    typedef Modifier :: Const<t_Object> :: t_constant t_constant;
    //! broacast instance
    BroadCast &_this;
    //! Constructor
    Object( BroadCast &__this ) : _this(__this) {}
    //! non constant serializing
    bool operator()( t_nonconstant &_ob )
      { ___OBJECTCODE }
    //! constant serializing
    bool operator()( t_constant &_ob )
      { ___OBJECTCODE }
  };
#endif
#endif

#undef ___TYPE__
#undef ___BUFF
#undef ___ENDBUFF
#undef ___BUFFNUM
#undef ___NUMTYPE
#undef ___CODE1
#undef ___CODE2
