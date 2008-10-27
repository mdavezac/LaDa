//
//  Version: $Id$
//
#ifndef  _DARWIN_ALLOYLAYERS_OBJECT_H_
#define  _DARWIN_ALLOYLAYERS_OBJECT_H_

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <string>

#include<opt/types.h>
#include<vff/layered.h>

#include "../bitstring.h"
#include "policies.h"

//! \cond
namespace Crystal
{
  class Structure;
}
namespace Vff
{
  class Keeper;
  types::t_real inplane_stress( const atat::rMatrix3d &_stress,
                                const atat::rVector3d &_dir );
}
//! \endcond

namespace GA
{
  namespace AlloyLayers
  {
    struct Object : public BitString::Object<>, 
                    public ::GA::Keepers::BandGap,
                    public ::GA::Keepers::OscStrength
    {
      //! The type of the BitString container
      typedef ::BitString::Object<> :: t_Container t_Container;
      //! Constructor
      Object() : ::BitString::Object<>(), ::GA::Keepers::BandGap(),
                 ::GA::Keepers::OscStrength() {}
      //! Copy Constructor
      Object   (const Object &_c)
             : ::BitString::Object<>(_c), ::GA::Keepers::BandGap(_c),
               ::GA::Keepers::OscStrength(_c) {};
      //! Loads from \a _node.
      bool Load( const TiXmlElement &_node )
        { return     ::GA::Keepers::BandGap::Load( _node ) 
                 and ::GA::Keepers::OscStrength::Load( _node ); }
      //! Saves to \a _node.
      bool Save( TiXmlElement &_node ) const
        { return     ::GA::Keepers::BandGap::Save( _node ) 
                 and ::GA::Keepers::OscStrength::Save( _node ); }
      //! Destructor
      virtual ~Object() {};
      //! Serializes a scalar individual.
      template<class Archive>
        void serialize(Archive & _ar, const unsigned int _version)
        {
          _ar & boost::serialization::base_object< ::BitString::Object<> >( *this ); 
          _ar & boost::serialization::base_object< ::GA::Keepers::BandGap >( *this ); 
          _ar & boost::serialization::base_object< ::GA::Keepers::OscStrength >( *this ); 
        }
    };

    //! \brief Old-style translation.
    //! \todo remove this function and replace it with translate.
    inline void operator<<( std::string &_str, const Object& _o )
      { Translate< Object > :: translate( _o, _str ); }
    //! \brief Old-style translation.
    //! \todo remove this function and replace it with translate.
    inline void operator<<( Object& _o, const std::string &_str )
      { Translate< Object > :: translate( _str, _o ); }
    //! \brief Old-style translation.
    //! \todo remove this function and replace it with translate.
    inline void operator<<( Crystal::Structure &_str, const Object& _o )
      { Translate< Object > :: translate( _o, _str ); }
    //! \brief Old-style translation.
    //! \todo remove this function and replace it with translate.
    inline void operator<<( Object& _o, const Crystal::Structure &_str )
      { Translate< Object > :: translate( _str, _o ); }

    //! Helper function which returns the inplane stress.
    template< class T_THIS >
     inline types::t_real inplane_stress( const Vff::Keeper& _o, const T_THIS &_this )
     {
       return Vff::inplane_stress( _o.stress, _this.get_direction() ); 
     }
  }
}

#endif
