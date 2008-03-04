//
//  Version: $Id: atom.h 548 2008-02-23 17:26:08Z davezac $
//
#ifndef _SERIALIZE_ATOM_H_


#ifndef ___TYPE__
#include <atat/serialize.h>

#define ___TYPE__ types::t_real
#include "atom.h"
#undef ___TYPE__

#define ___TYPE__ std::string
#include "atom.h"
#undef ___TYPE__

#define ___TYPE__ std::vector< std::string >
#include "atom.h"
#undef ___TYPE__
#define ___TYPE__ types::t_complex
#include "atom.h"
#undef ___TYPE__

#define _SERIALIZE_ATOM_H_

#else
  namespace mpi
  {
    //! Serializes atoms
    template <> struct BroadCast::Object< Ising_CE::Atom_Type<___TYPE__> >
    {
      //! Non constant type
      typedef Ising_CE::Atom_Type<___TYPE__> t_Object;
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
      { 
        return     _this.serialize( _ob.pos )
               and _this.serialize( _ob.type )
               and _this.serialize( _ob.freeze )
               and _this.serialize( _ob.site );
      }
      //! constant serializing
      bool operator()( t_constant &_ob )
      { 
        return     _this.serialize( _ob.pos )
               and _this.serialize( _ob.type )
               and _this.serialize( _ob.freeze )
               and _this.serialize( _ob.site );
      }
    };
  }
#endif
#endif
